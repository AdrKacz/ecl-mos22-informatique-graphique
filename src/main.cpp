#include <iostream>
#include <string>
#include <cmath>
#include <thread>
#include <random>
#include <omp.h>

#include "Camera/Camera.h"
#include "Vector/Vector.h"
#include "Ray/Ray.h"
#include "Sphere/Sphere.h"
#include "Mesh/TriangleMesh.h"
#include "Object/Object.h"
#include "Environment/Environment.h"
#include "Random/Random.h"


#define MESH_PATH "meshes/dog.obj"
#define USE_FOCAL_DISTANCE false
#define FOCAL_DISTANCE 20
#define USE_BOX_MULLER false
#define USE_INDIRECT_LIGHTING false
#define DYNAMIC_MOVEMENT false // true to move with arrow key, false to generate one image
#define INDIRECT_RAYS 1 // indirect rays number (no indirect ray with value 1)


#define HALF_BOX_DIMENSION 50 // in world unit
#define PIXELS_DIMENSION 128 // in pixels
#define OMP_NUM_THREADS 16

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../library/stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "../library/stb_image.h"

// ===== ===== ===== =====
// ===== ===== ===== ===== Helpers
// ===== ===== ===== =====

double clamp(const double& v, const double& low, const double& high) {
    return std::max(std::min(v, high), low);
}

// ===== ===== ===== =====
// ===== ===== ===== ===== Main
// ===== ===== ===== =====

void refresh(std::string filename, Environment& E, const int W, const int H, Camera& C, double I, double I_pow_factor, std::vector<unsigned char>& image) {
    auto start = std::chrono::high_resolution_clock::now();
    #pragma omp parallel for num_threads(OMP_NUM_THREADS) // schedule(dynamic, 1)
    for (int y = 0; y < H; y++) {
        for (int x = 0 ; x < W; x++) {
            Vector intensity = Vector();
            // INDIRECT_RAYS rays at random
            for (int k = 0; k < INDIRECT_RAYS; k++) {
                Vector u = Vector(x - W / 2 + 0.5, H - y - H / 2 + 0.5, C.z);
                if (USE_BOX_MULLER) {
                    Vector random_box = randh::box_muller();
                    u[0] += random_box[0];
                    u[1] += random_box[1];
                }
                u.normalize();

                Ray r;
                if (USE_FOCAL_DISTANCE) {
                    r = C.get_ray(u);
                } else {
                    r = Ray(C.get_position(), C.look_from(u));
                }
                Vector P, N;
                int object_index;
                if (E.intersect(r, P, N, &object_index))
                {
                    intensity = intensity + E.get_intensity(N, P, I, object_index, r, 0);
                }
            }
            intensity = intensity / INDIRECT_RAYS * 1.0;
            image[(y * W + x) * 3 + 0] = clamp(std::pow(intensity[0], I_pow_factor), 0., 255.);
            image[(y * W + x) * 3 + 1] = clamp(std::pow(intensity[1], I_pow_factor), 0., 255.);
            image[(y * W + x) * 3 + 2] = clamp(std::pow(intensity[2], I_pow_factor), 0., 255.);
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto diff_sec = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

    if (!DYNAMIC_MOVEMENT)
    {
        std::cout << "\nTemps pour la création de l'image: " << diff_sec.count() / 1000000.0 << "ms\n";
    }

    stbi_write_png(filename.c_str(), W, H, 3, &image[0], 0);
}

int main(int argc, char* argv[]) {
    // Get filename if any
    std::string filename = "image.png";
    if (argc > 1) {
        filename = std::string(argv[1]);
        if (filename.length() < 4 || filename.substr(filename.length() - 4) != ".png") {
            std::cout << "\e[1mError: Filename must end with \".png\"\e[0m\n";
            return 1;
        }
        filename = "outputs/" + filename;
    }

    // Dimension
	const int W = PIXELS_DIMENSION;
	const int H = PIXELS_DIMENSION;

    Environment E = Environment();

    // Objects
    std::cout << "Create mesh..." << std::endl;
    TriangleMesh mesh = TriangleMesh();
    std::cout << "Start reading mesh..." << std::endl;
    mesh.readOBJ(MESH_PATH);
    std::cout << "Create bounding box..." << std::endl;
    mesh.init_bounding_box();
    E.add_mesh(&mesh);
    // E.add_sphere(new Sphere(Vector(5, 10, 0), 2, Vector(1, 0, 0.5)));
    // E.add_sphere(new Sphere(Vector(0, 0, 0), 10, Vector(1, 0, 0.5)));
    // E.add_sphere(new Sphere(Vector(20, 0, 0), 5, Sphere::TYPE_REFLECTIVE));
    // E.add_sphere(new Sphere(Vector(5, 0, 20), 3, Sphere::TYPE_TRANSPARENT, 1.5));

    // Box
    const double wall_size = HALF_BOX_DIMENSION * 1e3;
    E.add_sphere(new Sphere(Vector(- wall_size - HALF_BOX_DIMENSION, 0, 0), wall_size)); // left
    E.add_sphere(new Sphere(Vector(+ wall_size + HALF_BOX_DIMENSION, 0, 0), wall_size)); // right
    E.add_sphere(new Sphere(Vector(0, - wall_size - HALF_BOX_DIMENSION, 0), wall_size)); // bottom
    E.add_sphere(new Sphere(Vector(0, + wall_size + HALF_BOX_DIMENSION, 0), wall_size)); // top
    E.add_sphere(new Sphere(Vector(0, 0, - wall_size - HALF_BOX_DIMENSION), wall_size)); // forward
    E.add_sphere(new Sphere(Vector(0, 0, + wall_size + HALF_BOX_DIMENSION), wall_size, Vector(0.8, 0.5, 0))); // backwards

    // Light
    // E.add_sphere(new Sphere(Vector(-30, 30, -30), 15, Sphere::TYPE_EMISSIVE, 1e9));
    // E.add_sphere(Sphere(Vector(30, 30, -30), 15, Sphere::TYPE_EMISSIVE, 1e9));
    E.add_light(Vector(30, 30, 25));
    // E.add_light(Vector(-20, -20, -15));

    // Indirect lighting
    E.use_indirect_lighting = USE_INDIRECT_LIGHTING;
    
    // Camera
    double alpha = 90 * M_PI / 180;
    Camera C = Camera(Vector(0, 0, 40), alpha, PIXELS_DIMENSION, FOCAL_DISTANCE);
	
	std::vector<unsigned char> image(W * H * 3, 0);

    double I = 1e10;
    double I_pow_factor = 1. / 2.2;

    std::string mvt_sequence = "";
    bool is_alive = true;
    char c = ' ';

    std::cout << "Start refreshing screen..." << std::endl;
    while (is_alive) {
        refresh(filename, E, W, H, C, I, I_pow_factor, image);
        if (!DYNAMIC_MOVEMENT)
        {
            is_alive = false;
            continue;
        }
        std::cin >> c;
        switch (c)
        {
        case 'z':
            C.move_forward(1.);
            break;
        case 's':
            C.move_forward(-1.);
            break;
        case 'd':
            C.move_right(1.);
            break;
        case 'q':
            C.move_right(-1.);
            break;
        case 'a':
            C.rotate(- M_PI / 12);
            break;
        case 'e':
            C.rotate(+ M_PI / 12);
            break;
        default:
            is_alive = false;
            break;
        }
        mvt_sequence += c;
    }
    
    if (DYNAMIC_MOVEMENT)
    {
        std::cout << "Movement Sequence:\n" << mvt_sequence; 
    }

	return 0;
}