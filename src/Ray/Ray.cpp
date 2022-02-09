#include "Ray.h"

Ray::Ray()
{
}

Ray::Ray(const Vector& a, const Vector& b)
{
    C = a;
    u = b;
}

Ray::~Ray()
{

}