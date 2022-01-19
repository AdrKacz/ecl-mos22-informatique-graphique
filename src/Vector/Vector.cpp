#include "Vector.h"

Vector::Vector()
{
    vec[0] = .0;
    vec[1] = .0;
    vec[2] = .0;
}

Vector::Vector(double x, double y, double z)
{
    vec[0] = x;
    vec[1] = y;
    vec[2] = z;
}

Vector::~Vector()
{
}

Vector Vector::operator+(const Vector& param) const
{
    return Vector(vec[0] + param[0], vec[1] + param[1], vec[2] + param[2]);

}

Vector Vector::operator-(const Vector& param) const
{
    return Vector(vec[0] - param[0], vec[1] - param[1], vec[2] - param[2]);

}

Vector Vector::operator*(const double param) const
{
    return Vector(vec[0] * param, vec[1] * param, vec[2] * param);

}

Vector Vector::operator/(const double param) const
{
    return Vector(vec[0] / param, vec[1] / param, vec[2] / param);

}

double Vector::dot(const Vector& param) const {
    return vec[0] * param[0] + vec[1] * param[1] + vec[2] * param[2];
}

Vector Vector::cross(const Vector& param) {
    return Vector(vec[1] * param[2] - vec[2] * param[1], vec[2] * param[0] - vec[0] * param[2], vec[0] * param[1] - vec[1] * param[0]);
}

double Vector::norm2() const {
    return vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
}

double Vector::norm() const {
    return sqrt(norm2());
}

void Vector::normalize() {
    double n = norm();
    vec[0] /= n;
    vec[1] /= n;
    vec[2] /= n;
}