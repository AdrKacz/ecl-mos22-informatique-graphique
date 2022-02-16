#include <iostream>
#include <string>
#include <cmath>
#include <thread>
#include <random>

// #include <ncurses.h>

#include "src/Camera/Camera.h"
#include "src/Vector/Vector.h"
#include "src/Ray/Ray.h"
#include "src/Sphere/Sphere.h"
#include "src/Mesh/TriangleMesh.h"
#include "src/Object/Object.h"
#include "src/Environment/Environment.h"


#define MESH_PATH "../meshes/dog.obj"
#define USE_FOCAL_DISTANCE false
#define FOCAL_DISTANCE 40
#define USE_BOX_MULLER true
#define BOX_MULLER_SIGMA 0.5
#define USE_INDIRECT_LIGHTING false
#define DYNAMIC_MOVEMENT false // true to move with arrow key, false to generate one image
#define INDIRECT_RAYS 1 // indirect rays number


#define HALF_BOX_DIMENSION 50 // in world unit
#define PIXELS_DIMENSION 128 // in pixels
#define NB_THREAD_GRID 6 // grid n * n > NB_THREAD = NB_THREAD_GRID * NB_THREAD_GRID

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

std::default_random_engine rng;
std::uniform_real_distribution<double> r_unif(0.0, 1.0);
constexpr double epsilon = std::numeric_limits<double>::epsilon();
constexpr double two_pi = 2 * M_PI;
Vector box_muller() {
    double u, v;
    do
    {
        u = r_unif(rng);
    } while (u <= epsilon);
    v = r_unif(rng);

    double mag = BOX_MULLER_SIGMA * sqrt(-2.0 * log(u));
    double x = mag * cos(2 * M_PI * v);
    double y = mag * sin(2 * M_PI * v);

    return Vector(x, y, 0);
}

// ===== ===== ===== =====
// ===== ===== ===== ===== Main
// ===== ===== ===== =====

void refresh(std::string filename, std::thread threads[NB_THREAD_GRID * NB_THREAD_GRID], const int step, Environment& E, const int W, const int H, Camera& C, double I, double I_pow_factor, std::vector<unsigned char>& image) {
    auto start = std::chrono::high_resolution_clock::now();
    for (int x = 0; x < NB_THREAD_GRID; x++)    
    {
        for (int y = 0; y < NB_THREAD_GRID; y++) {
            // uncomment to see black diagonal (thread not active)
            // if (x == y) { continue; }
            threads[NB_THREAD_GRID * x + y] = std::thread([x, y, &step, &E, &W, &H, &C, &I, &I_pow_factor, &image]()
            {
                for (int i = y * step; i < std::min(H, (y + 1) * step); i++) {
                    for (int j = x * step; j < std::min(W, (x + 1) * step); j++) {
                        Vector intensity = Vector();
                        // INDIRECT_RAYS rays at random
                        for (int k = 0; k < INDIRECT_RAYS; k++) {
                            Vector u = Vector(j - W / 2 + 0.5, H - i - H / 2 + 0.5, C.z);
                            if (USE_BOX_MULLER) {
                                Vector random_box = box_muller();
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
                        image[(i*W + j) * 3 + 0] = clamp(std::pow(intensity[0], I_pow_factor), 0., 255.);
                        image[(i*W + j) * 3 + 1] = clamp(std::pow(intensity[1], I_pow_factor), 0., 255.);
                        image[(i*W + j) * 3 + 2] = clamp(std::pow(intensity[2], I_pow_factor), 0., 255.);
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

    if (!DYNAMIC_MOVEMENT)
    {
        std::cout << "Temps pour la création de l'image: " << diff_sec.count() / 1000000.0 << "ms\n";
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
    }
    filename = "../outputs/" + filename;

    // TODO: Random has always the same seed
    // Vector test_box = box_muller();
    // std::cout << test_box[0] << "\n";
    // std::cout << test_box[1] << "\n";
    // std::cout << test_box[2] << "\n";

    // Dimension
	const int W = PIXELS_DIMENSION;
	const int H = PIXELS_DIMENSION;

    // Thread
    const int step = std::ceil(PIXELS_DIMENSION * 1.0 / NB_THREAD_GRID * 1.0);
    std::thread threads[NB_THREAD_GRID * NB_THREAD_GRID];

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
	
	std::vector<unsigned char> image(W*H * 3, 0);

    double I = 1e10;
    double I_pow_factor = 1. / 2.2;

    std::string mvt_sequence = "";
    bool is_alive = true;
    char c = ' ';

    std::cout << "Start refreshing screen..." << std::endl;
    while (is_alive) {
        refresh(filename, threads, step, E, W, H, C, I, I_pow_factor, image);
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