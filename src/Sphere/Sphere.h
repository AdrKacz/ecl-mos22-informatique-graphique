#ifndef DEF_SPHERE
#define DEF_SPHERE

#include <cmath>
#include "../Ray/Ray.h"
#include "../Vector/Vector.h"

class Sphere
{
private:
    Vector O;
    double R;
public:
    Sphere();
    Sphere(const Vector&, double);
    ~Sphere();

    bool intersect(const Ray&);
    bool intersect(const Ray&, Vector&, Vector&);
    bool intersect(const Ray&, Vector&, Vector&, double&);
};

#endif
