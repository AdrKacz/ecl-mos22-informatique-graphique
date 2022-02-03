#ifndef DEF_CAMERA
#define DEF_CAMERA

#include "../Vector/Vector.h"
#include "../Ray/Ray.h"

class Camera
{
private: 
    Vector position = Vector();
    double angle = .0;
public:
    Camera();
    Camera(const Vector&);
    Camera(const Vector&, double, double);

    double focus_distance = 40.;
    double alpha = 0.;
    double z = 0.;
    const Vector get_position();
    const Ray get_ray(const Vector&);

    void move_forward(double);
    void move_right(double);
    void rotate(double);
    Vector look_from(const Vector&);

    ~Camera();
};


#endif