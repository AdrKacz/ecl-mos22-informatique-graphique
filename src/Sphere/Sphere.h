#ifndef DEF_SPHERE
#define DEF_SPHERE

#include <cmath>
#include "../Object/Object.h"
#include "../Ray/Ray.h"
#include "../Vector/Vector.h"

class Sphere : public Object
{
public:
    Vector O = Vector(0, 0, 5);
    double R = 1.;

    Sphere();
    Sphere(const Vector&, double);
    Sphere(const Vector&, double, const Vector&);
    Sphere(const Vector&, double, int);
    Sphere(const Vector&, double, int, double);
    ~Sphere();

    virtual bool intersect(const Ray&);
    virtual bool intersect(const Ray&, Vector&, Vector&);
    virtual bool intersect(const Ray&, Vector&, Vector&, double&);
};

#endif
