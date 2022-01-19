#include "Sphere.h"

Sphere::Sphere()
{
}

Sphere::Sphere(const Vector& a, double b)
{
    O = a;
    R = b;
}

Sphere::Sphere(const Vector& a, double b, const Vector& c)
{
    O = a;
    R = b;
    albedo = c;
}

Sphere::Sphere(const Vector& a, double b, int type)
{
    O = a;
    R = b;
    if (type == Sphere::TYPE_REFLECTIVE)
    {
        is_reflective = true;
    } else if (type == Sphere::TYPE_TRANSPARENT) 
    {
        is_transparent = true;
    }
    
}

Sphere::Sphere(const Vector& a, double b, int type, double c)
{
    O = a;
    R = b;
    if (type == Sphere::TYPE_REFLECTIVE)
    {
        is_reflective = true;
    } else if (type == Sphere::TYPE_TRANSPARENT) 
    {
        is_transparent = true;
        n = c;
    }
}

Sphere::~Sphere()
{
}

bool Sphere::intersect(const Ray& r) {
    double a = 1;
    double b = 2 * r.u.dot(r.C - O);
    double c = (r.C - O).norm2() - R * R;
    double delta = b * b - 4 * a * c;
    return delta >= 0;
}

bool Sphere::intersect(const Ray& r, Vector& P, Vector& N) {
    double T;
    return intersect(r, P, N, T);
}

bool Sphere::intersect(const Ray& r, Vector& P, Vector& N, double& T) {
    double a = 1;
    double b = 2 * r.u.dot(r.C - O);
    double c = (r.C - O).norm2() - R * R;
    double delta = b * b - 4 * a * c;

    if (delta > 0)
    {
        double sqrt_delta = sqrt(delta);
        double t1 = (-b - sqrt_delta) / 2 * a;
        double t2 = (-b + sqrt_delta) / 2 * a;
        if (t1 >= 0) {    
            T = t1;
            P = r.C + r.u * t1;
            N = (P - O);
            N.normalize();
            return true;
        } else if (t2 >= 0) {
            T = t2;
            P = r.C + r.u * t2;
            N = (P - O);
            N.normalize();
            return true;
        } else {
            return false;
        }
    } else if (delta == 0)
    {
        double sqrt_delta = sqrt(delta);
        double t = (-b - sqrt_delta) / 2 * a;
        if (t >= 0)
        {
            T = t;
            P = r.C + r.u * t;
            N = (P - O);
            N.normalize();
            return true;
        }
        else
        {
            return false;
        }
    } else {
        return false;
    }
}