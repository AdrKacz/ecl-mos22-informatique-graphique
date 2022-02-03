#ifndef DEF_ENVIRONMENT
#define DEF_ENVIRONMENT

// #include <iostream>

#include <cmath>
#include <vector>
#include "../Vector/Vector.h"
#include "../Ray/Ray.h"
#include "../Sphere/Sphere.h"

class Environment
{
private:
    std::vector<Sphere> spheres;
    std::vector<Vector> lights;
    std::vector<int> sphere_lights;
public:
    static const int BOUNCES_MAX = 5;
    static constexpr double N_AMBIANT = 1.;
    bool use_indirect_lighting = true;

    Environment();
    ~Environment();

    void add_sphere(const Sphere&);
    bool intersect(const Ray&);
    bool intersect(const Ray&, Vector&, Vector&);
    bool intersect(const Ray&, Vector&, Vector&, int*);

    int get_random_sphere_lights_index();

    void add_light(const Vector&);
    double get_intensity(const Vector&, const Vector&, const double&, const double&);
    Vector get_intensity(const Vector&, const Vector&, const double&, int);
    Vector get_intensity(const Vector&, const Vector&, const double&, int, const Ray&, int);
};

#endif