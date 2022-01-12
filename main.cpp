#include <iostream>
#include <string>
#include <cmath>
#include <thread>

#include "src/Vector/Vector.h"
#include "src/Ray/Ray.h"
#include "src/Sphere/Sphere.h"
#include "src/Environment/Environment.h"

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
	int W = 512;
	int H = 512;

    // Thread
    const int step = std::ceil(512 / NB_THREAD_GRID);
    std::thread threads[NB_THREAD_GRID * NB_THREAD_GRID];

    Environment E = Environment();
    E.add_sphere(Sphere(Vector(0, 0, 0), 10));

    E.add_sphere(Sphere(Vector(-10050, 0, 0), 10000));
    E.add_sphere(Sphere(Vector(+10050, 0, 0), 10000));
    E.add_sphere(Sphere(Vector(0, -10050, 0), 10000));
    E.add_sphere(Sphere(Vector(0, +10050, 0), 10000));
    E.add_sphere(Sphere(Vector(0, 0, -10050), 10000));
    E.add_sphere(Sphere(Vector(0, 0, +10050), 10000));

    E.add_light(Vector(30, 30, 25));
    // E.add_light(Vector(-20, -20, -15));
    

    Vector C = Vector(0, 0, 40);
    double alpha = 90 * M_PI / 180;
    double z = - W / (2 * tan(alpha / 2));
	
	std::vector<unsigned char> image(W*H * 3, 0);

    double rho = M_PI;
    double I = 5e6;

    auto start = std::chrono::high_resolution_clock::now();
    for (int x = 0; x < NB_THREAD_GRID; x++)    
    {
        for (int y = 0; y < NB_THREAD_GRID; y++) {
            // uncomment to see black diagonal (thread not active)
            // if (x == y) { continue; }
            threads[NB_THREAD_GRID * x + y] = std::thread([x, y, &step, &z, &E, &W, &H, &C, &rho, &I, &image]()
            {
                for (int i = y * step; i < std::min(H, (y + 1) * step); i++) {
                    for (int j = x * step; j < std::min(W, (x + 1) * step); j++) {    
                        Vector u = Vector(j - W / 2 + 0.5, H - i - H / 2 + 0.5, z);
                        u.normalize();
                        Ray r = Ray(C, u);

                        Vector P, N;
                        if (E.intersect(r, P, N))
                        {
                            double intensity = E.get_intensity(N, P, I, rho);
                            intensity = clamp(intensity, 0., 255.);
                            
                            image[(i*W + j) * 3 + 0] = intensity;
                            image[(i*W + j) * 3 + 1] = intensity;
                            image[(i*W + j) * 3 + 2] = intensity;
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