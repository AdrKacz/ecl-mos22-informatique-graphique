#ifndef DEF_ENVIRONMENT
#define DEF_ENVIRONMENT

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
public:
    Environment();
    ~Environment();

    void add_sphere(const Sphere&);
    bool intersect(const Ray&);
    bool intersect(const Ray&, Vector&, Vector&);

    void add_light(const Vector&);
    double get_intensity(const Vector&, const Vector&, const double&, const double&);
};

#endif