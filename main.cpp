#include <iostream>
#include <string>
#include <cmath>
#include <thread>

#define NB_THREAD_GRID 4 // grid n * n > NB_THREAD = NB_THREAD_GRID * NB_THREAD_GRID

#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "library/stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "library/stb_image.h"

double clamp(const double& v, const double& low, const double& high) {
    return std::max(std::min(v, high), low);
}

// ===== ===== ===== ===== Vector

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
    Vector operator/(const double) const;

    double dot(const Vector&) const;
    Vector cross(const Vector&);

    double norm2() const;
    double norm() const;
    void normalize();
};

// ===== ===== ===== ===== Ray

class Ray
{
private:
public:
    Vector C, u;
    Ray(const Vector&, const Vector&);
    ~Ray();
};

// ===== ===== ===== ===== Sphere

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

// ===== ===== ===== ===== Environment

class Environment
{
private:
    std::vector<Sphere> spheres;
    std::vector<Vector> lights;
public:
    Environment();
    ~Environment();

    void add_sphere(const Sphere&);
    bool intersect(const Ray&);
    bool intersect(const Ray&, Vector&, Vector&);

    void add_light(const Vector&);
    double get_intensity(const Vector&, const Vector&, const double&, const double&);
};

// ===== ===== ===== =====
// ===== ===== ===== ===== Vector
// ===== ===== ===== =====

Vector::Vector()
{
    vec[0] = 5;
    vec[1] = 5;
    vec[2] = 5;
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

// ===== ===== ===== =====
// ===== ===== ===== ===== Ray
// ===== ===== ===== =====

Ray::Ray(const Vector& a, const Vector& b)
{
    C = a;
    u = b;
}

Ray::~Ray()
{

}

// ===== ===== ===== =====
// ===== ===== ===== ===== Sphere
// ===== ===== ===== =====

Sphere::Sphere()
{
    O = Vector(0, 0, 5);
    R = 1;
}

Sphere::Sphere(const Vector& a, double b)
{
    O = a;
    R = b;
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
            }
        }
    }

    return has_intersected;  
}

double Environment::get_intensity(const Vector& N, const Vector& P, const double& I, const double& rho) {
    double intensity = 0.;
    for (int i = 0; i < lights.size(); i++)
    {
        Vector L = lights[i];
        Vector l = (L - P);
        l.normalize();
        double local_intensity = I * l.dot(N) * rho / (4 * M_PI * M_PI * (L - P).norm2());
        intensity += local_intensity;
    }
    return clamp(intensity, 0., 255.);
}

// ===== ===== ===== =====
// ===== ===== ===== ===== Main
// ===== ===== ===== =====

int main(int argc, char* argv[]) {
    // Get filename if any
    std::string filename = "image";
    if (argc > 1) {
        filename = std::string(argv[1]);
        if (filename.length() < 4 || filename.substr(filename.length() - 4) != ".png") {
            std::cout << "\e[1mError: Filename must end with \".png\"\e[0m\n";
            return 1;
        }
    }
    filename = "outputs/" + filename;

    // Dimension
	int W = 512;
	int H = 512;

    // Thread
    const int step = std::ceil(512 / NB_THREAD_GRID);
    std::thread threads[NB_THREAD_GRID * NB_THREAD_GRID];

    Environment E = Environment();
    E.add_sphere(Sphere(Vector(0, 0, 0), 10));

    E.add_sphere(Sphere(Vector(-10050, 0, 0), 10000));
    E.add_sphere(Sphere(Vector(+10050, 0, 0), 10000));
    E.add_sphere(Sphere(Vector(0, -10050, 0), 10000));
    E.add_sphere(Sphere(Vector(0, +10050, 0), 10000));
    E.add_sphere(Sphere(Vector(0, 0, -10050), 10000));
    E.add_sphere(Sphere(Vector(0, 0, +10050), 10000));

    E.add_light(Vector(30, 30, 25));
    // E.add_light(Vector(-20, -20, -15));
    

    Vector C = Vector(0, 0, 40);
    double alpha = 90 * M_PI / 180;
    double z = - W / (2 * tan(alpha / 2));
	
	std::vector<unsigned char> image(W*H * 3, 0);

    double rho = M_PI;
    double I = 5e6;

    auto start = std::chrono::high_resolution_clock::now();
    for (int x = 0; x < NB_THREAD_GRID; x++)    
    {
        for (int y = 0; y < NB_THREAD_GRID; y++) {
            // uncomment to see black diagonal (thread not active)
            // if (x == y) { continue; }
            threads[NB_THREAD_GRID * x + y] = std::thread([x, y, &step, &z, &E, &W, &H, &C, &rho, &I, &image]()
            {
                for (int i = y * step; i < std::min(H, (y + 1) * step); i++) {
                    for (int j = x * step; j < std::min(W, (x + 1) * step); j++) {    
                        Vector u = Vector(j - W / 2 + 0.5, H - i - H / 2 + 0.5, z);
                        u.normalize();
                        Ray r = Ray(C, u);

                        Vector P, N;
                        if (E.intersect(r, P, N))
                        {
                            double intensity = E.get_intensity(N, P, I, rho);
                            intensity = clamp(intensity, 0., 255.);
                            
                            image[(i*W + j) * 3 + 0] = intensity;
                            image[(i*W + j) * 3 + 1] = intensity;
                            image[(i*W + j) * 3 + 2] = intensity;
                        }
                    }
                }
            });
        }
    }
    
    for (int x = 0; x < NB_THREAD_GRID; x++)    
    {
        for (int y = 0; y < NB_THREAD_GRID; y++) {
            // uncomment to see black diagonal (thread not active)
            // if (x == y) { continue; }
            threads[x * NB_THREAD_GRID + y].join();
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto diff_sec = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

    std::cout << "Temps pour la création de l'image: " << diff_sec.count() / 1000000.0 << "ms\n";

    stbi_write_png(filename.c_str(), W, H, 3, &image[0], 0);
	

	return 0;
}