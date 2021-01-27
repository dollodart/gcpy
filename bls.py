"""
gcpy: python utilities for chromatographs
Copyright (C) 2021 David Ollodart <LICENSE>
"""
import numpy as np
from scipy.interpolate import interp1d


def lininterp_baseline_subtract(y, x, left_ips, right_ips):
    """

    To subtract bleed, which is characterized as an elongated peak, omit
    points which are found by peak-finding algorithms and then spline
    the remaining points, and subtract the interpolated values.

    """

    keep = np.ones(len(x), dtype=bool)
    for left, right in zip(left_ips, right_ips):
        keep[int(left):int(right)] = False

    f = interp1d(x[keep], y[keep])
    return y - f(x)
