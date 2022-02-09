#include "Environment.h"

// ===== ===== ===== =====
// ===== ===== ===== ===== Helpers
// ===== ===== ===== =====

#include <random>

#include <iostream>
#include <string>
static void debug(const std::string& s) {
    std::string sp = s + '\n';
    std::cout << sp;
}

std::default_random_engine generator;
std::uniform_real_distribution<double> distrib(0.0, 1.0);

static Vector random_cos(const Vector& N) {
    double r1 = distrib(generator);
    double r2 = distrib(generator);
    double x = cos(2 * M_PI * r1) * sqrt(1 - r2);
    double y = sin(2 * M_PI * r1) * sqrt(1 - r2);
    double z = sqrt(r2);

    Vector T1, T2;
    if (N[0] < N[1] && N[0] < N[2]) // Nx smallest
    {
        T1 = Vector(0, -N[2], N[1]);
    } else if (N[0] < N[1] && N[0] < N[2]) // Ny smallest
    {
        T1 = Vector(-N[2], 0, N[0]);
    } else { // nz smallest
        T1 = Vector(-N[1], N[0], 0);
    }
    T1.normalize();
    T2 = T1.cross(N);  
    T2.normalize();

    Vector w_i = T1 * x + T2 * y + N * z;
    w_i.normalize();
    return w_i;
}

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

void Environment::add_sphere(Sphere* s)
{
    objects.push_back(s);
    if (s->type == Object::TYPE_EMISSIVE) {
        sphere_lights.push_back(objects.size() - 1);
    }
}

void Environment::add_mesh(TriangleMesh* m)
{
    objects.push_back(m);
}

void Environment::add_light(const Vector& l) {
    lights.push_back(l);
}

bool Environment::intersect(const Ray& r)
{
    for (int i = 0; i < objects.size(); i++)
    {
        if (objects[i]->intersect(r)) {
            return true;
        }
    }
    return false;
    
}

int Environment::get_random_sphere_lights_index()
{
    int index = rand() % sphere_lights.size();
    return sphere_lights[index];
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
    for (int i = 0; i < objects.size(); i++)
    {
        double t;
        Vector p, n;
        if (objects[i]->intersect(r, p, n, t)) {
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
    return intensity / M_PI;
}

Vector Environment::get_intensity(const Vector& N, const Vector& P, const double& I, int object_index) {
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
        double local_intensity = I * l.dot(N) / (4 * M_PI * M_PI * (L - P).norm2());
        intensity += local_intensity;        
    }
    return objects[object_index]->albedo * intensity / M_PI;
}

Vector Environment::get_intensity(const Vector& N, const Vector& P, const double& I, int object_index, const Ray& bounce_i, int bounces) {
    if (bounces >= Environment::BOUNCES_MAX) {
        return Vector();
    }
    // Gather surface that receive light
    if (objects[object_index]->type == Object::TYPE_REFLECTIVE)
    {
        Vector r = bounce_i.u - N * N.dot(bounce_i.u) * 2;
        r.normalize();

        Ray bounce_r = Ray(P + N * .01, r);
        Vector bounce_P, bounce_N;
        if (intersect(bounce_r, bounce_P, bounce_N, &object_index)) {
            return get_intensity(bounce_N, bounce_P, I, object_index, bounce_r, bounces + 1);
        }
        return Vector();
    } else if (objects[object_index]->type == Object::TYPE_TRANSPARENT)
    {
        double n_ratio;
        Vector local_N;
        if (bounce_i.u.dot(N) < 0) // enter sphere
        {
            n_ratio = Environment::N_AMBIANT / objects[object_index]->n;
            local_N = N;
        } else // leave sphere
        {
            n_ratio = objects[object_index]->n / Environment::N_AMBIANT;
            local_N = N * -1.;
        }
        
        Vector t_N = local_N * std::sqrt(1. - std::pow(n_ratio, 2) * (1. - std::pow(bounce_i.u.dot(local_N), 2))) * -1.;
        Vector t_T = (bounce_i.u - local_N * local_N.dot(bounce_i.u)) * n_ratio;
        Vector t = t_N + t_T;
        t.normalize();

        Ray bounce_t = Ray(P + t * .01, t);
        Vector bounce_P, bounce_N;
        if (intersect(bounce_t, bounce_P, bounce_N, &object_index)) {
            return get_intensity(bounce_N, bounce_P, I, object_index, bounce_t, bounces + 1);
        }
        return Vector();
    }


    // Sphere is diffuse
    // Measure light from surface
    // Direct lighting
    double intensity = .0;
    for (int i = 0; i < lights.size(); i++)
    {
        Vector L = lights[i];
        Vector l = (L - P);
        l.normalize();
        double l_dot_N = l.dot(N);
        if (l_dot_N <= 0) {
            continue;
        }
        
        // Check if there is intersection between P and L
        double distance_to_light = (L - P).norm2();
        Ray local_r = Ray(P + N * .01, l);
        Vector local_P, local_N;
        if (intersect(local_r, local_P, local_N)) {
            // Check if object closer than light
            double distance_to_object = (local_P - P).norm2();
            if (distance_to_object <= distance_to_light) {
            //     if (object_index == 0 && intensity == 0) {
            //      debug(std::string("Shadow at :\t" + std::to_string(object_index)));
            //  } 
                continue;
            }
        }
        // TODO: Check formula (4*PI*PI or 4*PI)
        double local_intensity = I * l_dot_N / (4 * M_PI * distance_to_light);
        intensity += local_intensity;     
    } 
    // if (object_index == 0 && intensity == 0) {
    //         debug(std::string("\tShadow on :\t" + std::to_string(object_index)));
    // } 
    Vector punctual_light_color = objects[object_index]->albedo * intensity;
    if (!use_indirect_lighting) // if no indirect lighting, emissive sphere won't work (no ray to bouce up to them)
    {
        return punctual_light_color / M_PI;
    }

    // Indirect lighting via emissive light (use trick to avoid noise)
    Vector emissive_light_color;
    if (sphere_lights.size() > 0) {
        int random_sphere_light_index = get_random_sphere_lights_index();
        if (random_sphere_light_index != object_index) { // Well, not useful if you're the light you want indirect lighting from...
            // debug(std::string("Get Light from " + std::to_string(random_sphere_light_index) + " to " + std::to_string(sphere_index)));
            Sphere* random_emissive_sphere = dynamic_cast<Sphere*>(objects[random_sphere_light_index]);
            Vector l = (P - random_emissive_sphere->O);
            l.normalize();  
            Vector w_random = random_cos(l);
            Vector random_P_on_light = random_emissive_sphere->O + w_random * random_emissive_sphere->R;

            Vector w_l = (random_P_on_light - P);
            w_l.normalize();

            double w_l_dot_N = w_l.dot(N);
            if (w_l_dot_N > 0) { // Continue only if emissive sphere is more or less in front
                double distance_to_light = (random_P_on_light - P).norm2(); 
                Ray local_r = Ray(P + N * .01, w_l);
                Vector local_P, local_N;
                if (intersect(local_r, local_P, local_N)) { // Should ALWAYS intersect (at least with the emissive sphere)!
                    double distance_to_object = (local_P - P).norm2();
                    if (distance_to_object > distance_to_light * .99) { // .99 to avoid counting self-colliding
                        // Add contribution from emissive light
                        emissive_light_color = objects[object_index]->albedo * random_emissive_sphere->albedo \
                            * random_emissive_sphere->intensity * w_l_dot_N * (w_random * -1.).dot(N) \
                            / (4 * M_PI * distance_to_light * l.dot(w_random));
                    }
                } else {
                    // Should never happen
                    debug(std::string("[ERROR] Collided with nothing..."));
                }
            }
        }
    }

    // Indirect lighting via all the rest (including emissive light by chance)
    Vector w_i = random_cos(N);
    Ray indirect_r = Ray(P + N * .01, w_i);
    Vector indirect_P, indirect_N;
    Vector indirect_intensity;
    int object_index_indirect;
    if (intersect(indirect_r, indirect_P, indirect_N, &object_index_indirect)) {
        indirect_intensity = get_intensity(indirect_N, indirect_P, I, object_index_indirect, indirect_r, bounces + 1);
    }

    Vector color = (punctual_light_color + emissive_light_color) / M_PI + objects[object_index]->albedo * indirect_intensity;
    if (objects[object_index]->type == Object::TYPE_EMISSIVE)
    {
        Sphere* emissive_sphere = dynamic_cast<Sphere*>(objects[object_index]);
        return emissive_sphere->albedo * emissive_sphere->intensity / (4 * M_PI * emissive_sphere->R * emissive_sphere->R) + color;
    } else
    {
        return color;
    }
}
