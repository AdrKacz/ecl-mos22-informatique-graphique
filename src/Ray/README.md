# Ray

```c++
class Ray
{
private:
public:
    Vector C, u;
    Ray(const Vector&, const Vector&);
    ~Ray();
};
```

```c++
Ray::Ray(const Vector& a, const Vector& b)
{
    C = a;
    u = b;
}

Ray::~Ray()
{

}
```