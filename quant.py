from gcpy.gcfactors import retention_tree


def integrate(y):
    """
    Assumes uniformly spaced y.
    Scale with interval before or after integrating.
    """

    # first order good enough for high resolution

    return (y[0] + y[-1])/2 + y[1:-1].sum()


def assign_peaks(x, x_center, x_left, x_right, min_with_window=False):

    species = []
    n = len(x_center)
    for peakn in range(n):
        center = x_center[peakn]
        lb, ub = int(x_left[peakn] // 1), int(x_right[peakn] // 1)
        pwidth = x[ub] - x[lb]

        match = sorted(retention_tree[x[center]])
        if match:
            sp = match[0].data
            if min_with_window:
                mn = abs(match[0].length() - pwidth) +\
                    abs(x[center] - (match[0].begin + match[0].end) / 2)
                for rt in match:
                    wwidth = rt.length()
                    delta = abs(wwidth - pwidth) +\
                        abs(x[center] - (rt.begin + rt.end) / 2)
                    if delta < mn:
                        mn = delta
                        sp = rt.data
        else:
            sp = 'unknown'

        species.append(sp)

    return species


def integrate_peaks(y, x_left, x_right):
    n = len(x_left)
    areas = []
    for peakn in range(n):
        lb, ub = int(x_left[peakn] // 1), int(x_right[peakn] // 1)
        area = integrate(y[lb:ub])
        areas.append(area)
    return areas
