"""

The data structures used here can be changed to be more efficient. For example, you can
construct a binary search tree for intervals rather than iterate through
all of them to determine the closest match, as is done here.

"""

import numpy as np

species = ['A',
           'B',
           'C']

sens_facts = [0.5, 1, 1.5]

rtimes = np.array([
    [0.35, 0.4],
    [0.7, 0.8],
    [0.85, 0.9]
])

if __name__ == '__main__':  # create test data
    from rampy_funcs import gaussian
    import numpy.random as r
    import numpy as np
    from encode import encode
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
