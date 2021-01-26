from curve_funcs import gaussian
import numpy.random as r
import numpy as np
from encode import encode

from sensitivity_data import rtimes

r.seed(232015092)
xmin = 0
xmax = 1
n = 10000
x = np.linspace(xmin, xmax, n)
y = 0
for rt in rtimes:
    amplitude = r.random_sample() + 1.
    mean = rt.mean()
    std = r.random_sample() / 20.
    y += gaussian(x, amplitude, mean, std)
y += r.random(n) / 100
encode('data/test-data.CH', y)

y += r.random(n) / 10
encode('data/noisy-test-data.CH', y)

y += np.hstack( (np.linspace(0, 0.5, n // 2) **2, np.zeros(n // 2 + n % 2)) )
encode('data/baseline-and-noisy-test-data.CH', y)
