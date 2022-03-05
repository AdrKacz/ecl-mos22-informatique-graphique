#ifndef RANDOMHELPER
#define RANDOMHELPER

#include <random>
#include "../Vector/Vector.h"

#define BOX_MULLER_SIGMA 0.5

#ifndef NUM_THREADS
#define NUM_THREADS 16
#endif

namespace randh
{
    Vector box_muller();
    Vector random_cos(const Vector &N);
}

#endif