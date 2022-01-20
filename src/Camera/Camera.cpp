#include "Camera.h"

Camera::Camera()
{
}

Camera::Camera(const Vector& p)
{
    position = p;
}

Camera::Camera(const Vector& p, double a)
{
    position = p;
    angle = a;
}

Camera::~Camera()
{
}

const Vector Camera::get_position()
{
    return position;
}

void Camera::move_forward(double d)
{
    Vector direction_xz = Vector(0, 0, -1);
    direction_xz.normalize();
    position = position + look_from(direction_xz) * 2. * d;
}

void Camera::rotate(double a)
{
    angle += a;
    while (angle > 2 * M_PI) {
        angle -= 2 * M_PI;
    }

    while (angle < 0) {
        angle += 2 * M_PI;
    }
}

Vector Camera::look_from(const Vector& v)
{
    double x = v[0];
    double z = v[2];
    double cos_a = std::cos(angle);
    double sin_a = std::sin(angle);

    Vector v_oriented = Vector(x * cos_a - z * sin_a, v[1], x * sin_a + z * cos_a);
    v_oriented.normalize();

    return v_oriented;
}