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
