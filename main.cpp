#include <iostream>
#include <string>
#include <cmath>
#include <thread>

#include "src/Vector/Vector.h"
#include "src/Ray/Ray.h"
#include "src/Sphere/Sphere.h"
#include "src/Environment/Environment.h"

#define HALF_BOX_DIMENSION 50 // in world unit
#define PIXELS_DIMENSION 256 // in pixels
#define NB_THREAD_GRID 4 // grid n * n > NB_THREAD = NB_THREAD_GRID * NB_THREAD_GRID

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "library/stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "library/stb_image.h"

// ===== ===== ===== =====
// ===== ===== ===== ===== Helpers
// ===== ===== ===== =====

double clamp(const double& v, const double& low, const double& high) {
    return std::max(std::min(v, high), low);
}

// ===== ===== ===== =====
// ===== ===== ===== ===== Main
// ===== ===== ===== =====

int main(int argc, char* argv[]) {
    // Get filename if any
    std::string filename = "image.png";
    if (argc > 1) {
        filename = std::string(argv[1]);
        if (filename.length() < 4 || filename.substr(filename.length() - 4) != ".png") {
            std::cout << "\e[1mError: Filename must end with \".png\"\e[0m\n";
            return 1;
        }
    }
    filename = "../outputs/" + filename;

    // Dimension
	int W = PIXELS_DIMENSION;
	int H = PIXELS_DIMENSION;

    // Thread
    const int step = std::ceil(PIXELS_DIMENSION / NB_THREAD_GRID);
    std::thread threads[NB_THREAD_GRID * NB_THREAD_GRID];

    Environment E = Environment();
    E.add_sphere(Sphere(Vector(0, 0, 0), 10, Vector(1, 0, 0.5)));
    E.add_sphere(Sphere(Vector(20, 0, 0), 5, Sphere::TYPE_REFLECTIVE));
    E.add_sphere(Sphere(Vector(5, 0, 20), 3, Sphere::TYPE_TRANSPARENT, 1.5));

    // Box
    const double wall_size = HALF_BOX_DIMENSION * 1e3;
    E.add_sphere(Sphere(Vector(- wall_size - HALF_BOX_DIMENSION, 0, 0), wall_size));
    E.add_sphere(Sphere(Vector(+ wall_size + HALF_BOX_DIMENSION, 0, 0), wall_size));
    E.add_sphere(Sphere(Vector(0, - wall_size - HALF_BOX_DIMENSION, 0), wall_size));
    E.add_sphere(Sphere(Vector(0, + wall_size + HALF_BOX_DIMENSION, 0), wall_size));
    E.add_sphere(Sphere(Vector(0, 0, - wall_size - HALF_BOX_DIMENSION), wall_size));
    E.add_sphere(Sphere(Vector(0, 0, + wall_size + HALF_BOX_DIMENSION), wall_size, Vector(0.8, 0.5, 0)));

    // Light
    E.add_light(Vector(30, 30, 25));
    // E.add_light(Vector(-20, -20, -15));
    

    Vector C = Vector(0, 0, 40);
    double alpha = 90 * M_PI / 180;
    double z = - W / (2 * tan(alpha / 2));
	
	std::vector<unsigned char> image(W*H * 3, 0);

    double rho = M_PI;
    double I = 1e10;
    double I_pow_factor = 1. / 2.2;

    auto start = std::chrono::high_resolution_clock::now();
    for (int x = 0; x < NB_THREAD_GRID; x++)    
    {
        for (int y = 0; y < NB_THREAD_GRID; y++) {
            // uncomment to see black diagonal (thread not active)
            // if (x == y) { continue; }
            threads[NB_THREAD_GRID * x + y] = std::thread([x, y, &step, &z, &E, &W, &H, &C, &rho, &I, &I_pow_factor, &image]()
            {
                for (int i = y * step; i < std::min(H, (y + 1) * step); i++) {
                    for (int j = x * step; j < std::min(W, (x + 1) * step); j++) {    
                        Vector u = Vector(j - W / 2 + 0.5, H - i - H / 2 + 0.5, z);
                        u.normalize();
                        Ray r = Ray(C, u);

                        Vector P, N;
                        int sphere_index;
                        if (E.intersect(r, P, N, &sphere_index))
                        {
                            Vector intensity = E.get_intensity(N, P, I, rho, sphere_index, r, 0);
                            image[(i*W + j) * 3 + 0] = clamp(std::pow(intensity[0], I_pow_factor), 0., 255.);;
                            image[(i*W + j) * 3 + 1] = clamp(std::pow(intensity[1], I_pow_factor), 0., 255.);;
                            image[(i*W + j) * 3 + 2] = clamp(std::pow(intensity[2], I_pow_factor), 0., 255.);;
                        }
                    }
                }
            });
        }
    }
    
    for (int x = 0; x < NB_THREAD_GRID; x++)    
    {
        for (int y = 0; y < NB_THREAD_GRID; y++) {
            // uncomment to see black diagonal (thread not active)
            // if (x == y) { continue; }
            threads[x * NB_THREAD_GRID + y].join();
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto diff_sec = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

    std::cout << "Temps pour la création de l'image: " << diff_sec.count() / 1000000.0 << "ms\n";

    stbi_write_png(filename.c_str(), W, H, 3, &image[0], 0);
	

	return 0;
}