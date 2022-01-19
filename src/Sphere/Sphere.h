#ifndef DEF_SPHERE
#define DEF_SPHERE

#include <cmath>
#include "../Ray/Ray.h"
#include "../Vector/Vector.h"

class Sphere
{
private:
    Vector O = Vector(0, 0, 5);
    double R = 1.;
public:
    static const int TYPE_DIFFUSE = 0;
    static const int TYPE_REFLECTIVE = 1;
    static const int TYPE_TRANSPARENT = 2;

    bool is_reflective = false;
    Vector albedo = Vector(1., 1., 1.);
    bool is_transparent = false;
    double n = 1.;
    
    Sphere();
    Sphere(const Vector&, double);
    Sphere(const Vector&, double, const Vector&);
    Sphere(const Vector&, double, int);
    Sphere(const Vector&, double, int, double);
    ~Sphere();

    bool intersect(const Ray&);
    bool intersect(const Ray&, Vector&, Vector&);
    bool intersect(const Ray&, Vector&, Vector&, double&);
};

#endif
