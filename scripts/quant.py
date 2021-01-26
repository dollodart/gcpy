from gcpy.decode import read_data
from gcpy.bls import lininterp_baseline_subtract
from gcpy.gcfactors import retention_tree, sensitivity_factors
from scipy.signal import find_peaks, peak_widths
import matplotlib.pyplot as plt
import numpy as np

# read in data
data = read_data('data/test-data.CH')
x = data['time']
y = data['tic']

# this subtracts the baseline by linearly interpolating between all points
# not identified as being part of peaks
peak_centers, peak_prop = find_peaks(y, width=40, height=1, rel_height=0.99)
y2 = lininterp_baseline_subtract(
    y, x, peak_prop['left_ips'], peak_prop['right_ips'])

# plot to visually inspect baseline subtraction
plt.plot(x, y, '--', label='orig')
plt.plot(x, y2, label='orig - bls')
for i in range(len(peak_centers)):
    plt.plot([x[int(peak_prop['left_ips'][i] // 1)],
              x[int(peak_prop['right_ips'][i] // 1)]], [y[peak_centers[i]]] * 2, 'b-|')
plt.legend()
plt.show()

y = y2

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

areas = integrate_peaks(y, peak_prop['left_ips'], peak_prop['right_ips'])
ids = assign_peaks(x, peak_centers, peak_prop['left_ips'], peak_prop['right_ips'])

wareas = [areas[i] / sensitivity_factors[ids[i]] for i in range(len(areas))]
total_warea = sum(wareas)

print('Species\tcomposition(mol%)')
for i in range(len(ids)):
    comp = 100 * wareas[i] / total_warea
    if comp > 0.1:
        print(f'{ids[i]}\t{comp:.2f}')
