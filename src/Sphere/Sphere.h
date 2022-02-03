#ifndef DEF_SPHERE
#define DEF_SPHERE

#include <cmath>
#include "../Ray/Ray.h"
#include "../Vector/Vector.h"

class Sphere
{
public:
    static const int TYPE_DIFFUSE = 0;
    static const int TYPE_REFLECTIVE = 1;
    static const int TYPE_TRANSPARENT = 2;
    static const int TYPE_EMISSIVE = 3;

    Vector O = Vector(0, 0, 5);
    double R = 1.;

    int type = 0.; // default to diffuse

    Vector albedo = Vector(1., 1., 1.); // diffuse and emissive
    double n = 1.; // transparent
    double intensity = 1.; // emissive
    
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
