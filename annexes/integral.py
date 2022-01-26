from math import cos, sin, pow, exp, sqrt, log, pi
from random import uniform

sigma = 0.5
mu = 0

f = lambda x: pow(cos(x), 10) if -pi / 2 <= x <= pi / 2 else 0
p = lambda x : exp(- pow(x - mu, 2) / pow(sigma, 2) * .5) / (sigma * sqrt(2 * pi))

box_muller = lambda u, v : {
    'x': cos(2 * pi * u) * sigma * sqrt(-2 * log(v)) + mu,
    'y': sin(2 * pi * u) * sigma * sqrt(-2 * log(v)) + mu,
}

# 1D
N = int(1e6)
I = sum([f(x) / p(x) for x in map(lambda _ : box_muller(uniform(0, 1), uniform(0, 1))['x'], range(N))]) / N
print(I) # 0.7727555623750062

# 3D
N = int(1e6)
I = sum([f(x) * f(y) * f(z) / p(x) / p(y) / p(z) for \
    x, y, z in map(lambda _ : ( \
        box_muller(uniform(0, 1), uniform(0, 1))['x'], \
        box_muller(uniform(0, 1), uniform(0, 1))['x'], \
        box_muller(uniform(0, 1), uniform(0, 1))['x']), \
    range(N))]) / N
print(I) # 0.46155941321408395