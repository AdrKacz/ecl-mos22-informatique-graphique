#include "Environment.h"

#include <iostream>

// ===== ===== ===== =====
// ===== ===== ===== ===== Helpers
// ===== ===== ===== =====

static double clamp(const double& v, const double& low, const double& high) {
    return std::max(std::min(v, high), low);
}

// ===== ===== ===== =====
// ===== ===== ===== ===== Environment
// ===== ===== ===== =====

Environment::Environment()
{
}

Environment::~Environment()
{
}

void Environment::add_sphere(const Sphere& s)
{
    spheres.push_back(s);
}

void Environment::add_light(const Vector& l) {
    lights.push_back(l);
}

bool Environment::intersect(const Ray& r)
{
    for (int i = 0; i < spheres.size(); i++)
    {
        if (spheres[i].intersect(r)) {
            return true;
        }
    }
    return false;
    
}

bool Environment::intersect(const Ray& r, Vector& P, Vector& N)
{
    int i;
    return intersect(r, P, N, &i);  
}

bool Environment::intersect(const Ray& r, Vector& P, Vector& N, int* index)
{
    bool has_intersected = false;
    double T = std::numeric_limits<double>::max();
    for (int i = 0; i < spheres.size(); i++)
    {
        double t;
        Vector p, n;
        if (spheres[i].intersect(r, p, n, t)) {
            has_intersected = true;
            if (t < T) {
                T = t;
                P = p;
                N = n;
                *index = i;
            }
        }
    }

    return has_intersected; 
}

double Environment::get_intensity(const Vector& N, const Vector& P, const double& I, const double& rho) {
    double intensity = .0;
    for (int i = 0; i < lights.size(); i++)
    {
        Vector L = lights[i];
        Vector l = (L - P);
        l.normalize();
        // Check if there is intersection between P and L
        Ray local_r = Ray(P + N * .01, l);
        Vector local_P, local_N;
        if (intersect(local_r, local_P, local_N)) {
            // Check if object closer than light
            double distance_to_object = (local_P - P).norm2();
            double distance_to_light = (L - P).norm2();
            if (distance_to_object <= distance_to_light) {
                continue;
            }
        }
        double local_intensity = I * l.dot(N) * rho / (4 * M_PI * M_PI * (L - P).norm2());
        intensity += local_intensity;
        
        
    }
    return intensity;
}

Vector Environment::get_intensity(const Vector& N, const Vector& P, const double& I, const double& rho, int sphere_index) {
    Vector intensity = Vector();
    for (int i = 0; i < lights.size(); i++)
    {
        Vector L = lights[i];
        Vector l = (L - P);
        l.normalize();
        // Check if there is intersection between P and L
        Ray local_r = Ray(P + N * .01, l);
        Vector local_P, local_N;
        if (intersect(local_r, local_P, local_N)) {
            // Check if object closer than light
            double distance_to_object = (local_P - P).norm2();
            double distance_to_light = (L - P).norm2();
            if (distance_to_object <= distance_to_light) {
                continue;
            }
        }
        double local_intensity = I * l.dot(N) * rho / (4 * M_PI * M_PI * (L - P).norm2());
        intensity = intensity + spheres[sphere_index].albedo * local_intensity;        
    }
    return intensity;
}

Vector Environment::get_intensity(const Vector& N, const Vector& P, const double& I, const double& rho, int sphere_index, const Ray& bounce_i, int bounces) {
    if (bounces >= Environment::BOUNCES_MAX) {
        return Vector();
    }
    
    // Gather surface that receive light
    if (spheres[sphere_index].is_reflective)
    {
        Vector r = bounce_i.u - N * N.dot(bounce_i.u) * 2;
        r.normalize();

        Ray bounce_r = Ray(P + N * .01, r);
        Vector bounce_P, bounce_N;
        if (intersect(bounce_r, bounce_P, bounce_N, &sphere_index)) {
            return get_intensity(bounce_N, bounce_P, I, rho, sphere_index, bounce_r, bounces + 1);
        }
        return Vector();
    } else if (spheres[sphere_index].is_transparent)
    {
        double n_ratio;
        Vector local_N;
        if (bounce_i.u.dot(N) < 0) // enter sphere
        {
            n_ratio = Environment::N_AMBIANT / spheres[sphere_index].n;
            local_N = N;
        } else // leave sphere
        {
            n_ratio = spheres[sphere_index].n / Environment::N_AMBIANT;
            local_N = N * -1.;
        }
        
        Vector t_N = local_N * std::sqrt(1. - std::pow(n_ratio, 2) * (1. - std::pow(bounce_i.u.dot(local_N), 2))) * -1.;
        Vector t_T = (bounce_i.u - local_N * local_N.dot(bounce_i.u)) * n_ratio;
        Vector t = t_N + t_T;
        t.normalize();

        Ray bounce_t = Ray(P + t * .01, t);
        Vector bounce_P, bounce_N;
        if (intersect(bounce_t, bounce_P, bounce_N, &sphere_index)) {
            return get_intensity(bounce_N, bounce_P, I, rho, sphere_index, bounce_t, bounces + 1);
        }
        return Vector();
    }

    // Measure light from surface
    Vector intensity = Vector();
    for (int i = 0; i < lights.size(); i++)
    {
        Vector L = lights[i];
        Vector l = (L - P);
        l.normalize();
        // Check if there is intersection between P and L
        Ray local_r = Ray(P + N * .01, l);
        Vector local_P, local_N;
        if (intersect(local_r, local_P, local_N)) {
            // Check if object closer than light
            double distance_to_object = (local_P - P).norm2();
            double distance_to_light = (L - P).norm2();
            if (distance_to_object <= distance_to_light) {
                continue;
            }
        }
        double local_intensity = I * l.dot(N) * rho / (4 * M_PI * M_PI * (L - P).norm2());
        intensity = intensity + spheres[sphere_index].albedo * local_intensity;        
    }
    return intensity;
}