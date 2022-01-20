#ifndef DEF_CAMERA
#define DEF_CAMERA

#include "../Vector/Vector.h"

class Camera
{
private: 
    Vector position = Vector();
    double angle = .0;
public:
    Camera();
    Camera(const Vector&);
    Camera(const Vector&, double);

    const Vector get_position();

    void move_forward(double);
    void rotate(double);
    Vector look_from(const Vector&);

    ~Camera();
};


#endif