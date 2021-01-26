import numpy as np
from scipy.interpolate import interp1d

def lininterp_baseline_subtract(y, x, left_ips, right_ips):
    """

    To subtract bleed, which is characterized as an elongated peak, omit
    points which are found by peak-finding algorithms and then spline
    the remaining points, and subtract the interpolated values.

    """

    keep = np.ones(len(x), dtype=bool)
    for l, r in zip(left_ips, right_ips):
        l = int(l)
        r = int(r)
        keep[l:r] = False

    f = interp1d(x[keep], y[keep])
    return y - f(x)


if __name__ == '__main__':

    from py_chemplexity import read_data
    import matplotlib.pyplot as plt

    # read in data
    data = read_data('data/baseline-and-noisy-test-data.CH')
    x = data['time']
    y = data['tic']
    als_smoothing_inspection(y, l_fineness=5, p_fineness=3, log10l_lims=[1, 5], log10p_lims=[-1, -0.301])
    plt.show()
