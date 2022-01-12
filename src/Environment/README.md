# Environment

```c++
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
}
```

```c++
Environment::Environment()
{
}

Environment::~Environment()
{
}
```

```c++
void Environment::add_sphere(const Sphere& s)
{
    spheres.push_back(s);
}

void Environment::add_light(const Vector& l) {
    lights.push_back(l);
}
```

```c++
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
```

```c++
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
```