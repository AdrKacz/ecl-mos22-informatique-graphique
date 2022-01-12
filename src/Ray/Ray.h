#ifndef DEF_RAY
#define DEF_RAY

#include "../Vector/Vector.h"

class Ray
{
private:
public:
    Vector C, u;
    Ray(const Vector&, const Vector&);
    ~Ray();
};

#endif