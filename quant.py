from decode import read_data
from bls import lininterp_baseline_subtract
from scipy.signal import find_peaks, peak_widths
from scipy.integrate import trapz
import matplotlib.pyplot as plt
import numpy as np
from sensitivity_data import species, rtimes, sens_facts

# read in data
data = read_data('data/test-data.CH')
x = data['time']
y = data['tic']

# this subtracts the baseline by linearly interpolating between all points
# not identified as being part of peaks
peak_pos, peak_prop = find_peaks(y, width=40, height=1, rel_height=0.99)
y2 = lininterp_baseline_subtract(
    y, x, peak_prop['left_ips'], peak_prop['right_ips'])

# plot to visually inspect baseline subtraction
plt.plot(x, y, '--', label='orig')
plt.plot(x, y2, label='orig - bls')
for i in range(len(peak_pos)):
    plt.plot([x[int(peak_prop['left_ips'][i] // 1)],
              x[int(peak_prop['right_ips'][i] // 1)]], [y[peak_pos[i]]] * 2, 'b-|')
plt.legend()
plt.show()

y = y2

peaks = zip(peak_pos, peak_prop['left_ips'], peak_prop['right_ips'])
peaks_lst = []
total = 0
for bound in peaks:
    center, lb, ub = bound
    lb, ub = int(lb // 1), int(ub // 1)
    integral = trapz(y[lb:ub])
    time = x[center]
    diffs = []
    for index, window in enumerate(rtimes):
        if time - window[0] > 0 and window[1] - time > 0:
            diffs.append([abs(time - window.mean()), index])
    if len(diffs) > 0:
        best_match = min(diffs, key=lambda x: x[0])
        index = best_match[1]
        sp = species[index]
        sf = sens_facts[index]
    else:
        sp = 'unknown'
        sf = 1

    peaks_lst.append([sp, integral / sf])
    total += integral / sf

print('Species\tcomposition(mol%)')
for species, area in peaks_lst:
    comp = 100 * area / total
    if comp > 0.1:  # greater than 1 part in 1000
        print(f'{species}\t{comp:.2f}')
