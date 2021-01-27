"""
gcpy: python utilities for chromatographs
Copyright (C) 2021 David Ollodart <LICENSE>
"""

from intervaltree import Interval, IntervalTree
retention_tree = IntervalTree()
sensitivity_factors = dict(unknown=1)

for lb, ub, sp, sf in ((0.35, 0.50, 'A', 0.5),
                       (0.70, 0.80, 'B', 1.0),
                       (0.85, 0.90, 'C', 1.5)):
    retention_tree.add(Interval(lb, ub, data=sp))
    sensitivity_factors[sp] = sf
