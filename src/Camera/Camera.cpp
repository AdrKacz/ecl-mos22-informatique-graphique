#include "Camera.h"

// ===== ===== ===== =====
// ===== ===== ===== ===== Helpers
// ===== ===== ===== =====

#include <random>

std::default_random_engine rng_camera;
std::uniform_real_distribution<double> r_unif_camera(0.0, 1.0);

// ===== ===== ===== =====
// ===== ===== ===== ===== Camera
// ===== ===== ===== =====

Camera::Camera()
{
}

Camera::Camera(const Vector& p)
{
    position = p;
}

Camera::Camera(const Vector& p, double a, double pixel_dim)
{
    position = p;
    alpha = a;
    z = - pixel_dim / (2. * tan(alpha / 2.));
}

Camera::Camera(const Vector& p, double a, double pixel_dim, double focal)
{
    position = p;
    alpha = a;
    z = - pixel_dim / (2. * tan(alpha / 2.));
    focus_distance = focal;
}

Camera::~Camera()
{
}

const Vector Camera::get_position()
{
    return position;
}

const Ray Camera::get_ray(const Vector& u)
{
    // TODO: Circular aperture, Ray r = Ray(C.get_position(), C.look_from(u));
    double dx = (r_unif_camera(rng_camera) - .5) * 5.;
    double dy = (r_unif_camera(rng_camera) - .5) * 5.;

    Vector origin = position + Vector(dx, dy, 0.);
    Vector point_to = position + look_from(u) * focus_distance - origin;
    point_to.normalize();

    return Ray(origin, point_to);
}

void Camera::move_forward(double d)
{
    Vector direction_xz = Vector(0, 0, -1);
    direction_xz.normalize();
    position = position + look_from(direction_xz) * 2. * d;
}

void Camera::move_right(double d)
{
    Vector direction_xz = Vector(1, 0, 0);
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