from gcpy.gcfactors import retention_tree


def integrate(y):
    """
    Assumes uniformly spaced y--scale with interval before integrating.
    """

    # first order good enough for high resolution

    return (y[0] + y[-1])/2 + y[1:-1].sum() 


def assign_peaks(x, x_center, x_left, x_right):

    species = []
    n = len(x_center)
    for peakn in range(n):
        center = x_center[peakn]
        lb, ub = int(x_left[peakn] // 1), int(x_right[peakn] // 1)
        t = x[center]

        match = sorted(retention_tree[x[center]])
        if match:
            sp = match[0].data
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
