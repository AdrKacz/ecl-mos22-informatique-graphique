#ifndef DEF_VECTOR
#define DEF_VECTOR

#include <cmath>
#include <string>

class Vector
{
private:
    double vec[3];
public:
    Vector();
    Vector(double x, double y, double z);
    ~Vector();

    double operator[](int i) const { return vec[i]; };
    double& operator[](int i) {return vec[i];};

    Vector operator+(const Vector&) const;
    Vector operator-(const Vector&) const;
    Vector operator*(const double) const;
    Vector operator*(const Vector param) const;
    Vector operator/(const double) const;

    double dot(const Vector&) const;
    Vector cross(const Vector&);

    double norm2() const;
    double norm() const;
    void normalize();

    std::string to_string();
};

#endif