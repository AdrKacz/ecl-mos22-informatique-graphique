#include <random>
#include <iostream>
#include <omp.h>
#include <cmath>
#include "../Vector/Vector.h"
#include "Random.h"

std::default_random_engine engine[NUM_THREADS];
std::uniform_real_distribution<double> uniform(0.0, 1.0);

constexpr double epsilon = std::numeric_limits<double>::epsilon();
constexpr double two_pi = 2 * M_PI;

Vector randh::box_muller()
{
    int thread_num = omp_get_thread_num();

    double u, v;
    do
    {
        u = uniform(engine[thread_num]);
    } while (u <= epsilon);
    v = uniform(engine[thread_num]);

    double mag = BOX_MULLER_SIGMA * sqrt(-2.0 * log(u));
    double x = mag * cos(2 * M_PI * v);
    double y = mag * sin(2 * M_PI * v);

    return Vector(x, y, 0);
};

Vector randh::random_cos(const Vector &N) // Retourne omega_i
{
    int thread_num = omp_get_thread_num();

    double r1 = uniform(engine[thread_num]);
    double r2 = uniform(engine[thread_num]);
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