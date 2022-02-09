#ifndef DEF_OBJECT
#define DEF_OBJECT

#include "../Vector/Vector.h"
#include "../Ray/Ray.h"

class Object
{
public:
    int type = 0.; // default to diffuse

    Vector albedo = Vector(1., 1., 1.); // diffuse and emissive
    double n = 1.; // transparent
    double intensity = 1.; // emissive

    static const int TYPE_DIFFUSE = 0;
    static const int TYPE_REFLECTIVE = 1;
    static const int TYPE_TRANSPARENT = 2;
    static const int TYPE_EMISSIVE = 3;

    virtual bool intersect(const Ray&) = 0;
    virtual bool intersect(const Ray&, Vector&, Vector&) = 0;
    virtual bool intersect(const Ray&, Vector&, Vector&, double&) = 0;

    Object();
    ~Object();
};

#endif