from gcpy.decode import read_data
from gcpy.bls import lininterp_baseline_subtract
from gcpy.gcfactors import sensitivity_factors
from gcpy.quant import integrate_peaks, assign_peaks
from scipy.signal import find_peaks
import matplotlib.pyplot as plt

# read in data
data = read_data('data/test-data.CH')
#data = read_data('data/noisy-test-data.CH')
#data = read_data('data/baseline-and-noisy-test-data.CH')
x = data['time']
y = data['tic']
#from gcpy.smooth import whittaker_smooth
#y = whittaker_smooth(y, 1e3)

# this subtracts the baseline by linearly interpolating between all points
# not identified as being part of peaks
peak_centers, peak_prop = find_peaks(y, width=40, height=1, rel_height=0.99)
peak_lefts, peak_rights = peak_prop['left_ips'], peak_prop['right_ips']

# for noisy and baseline data
#from scipy.signal import peak_widths
#peak_centers, peak_prop = find_peaks(y, height=1,
#                                    prominence=0.5,rel_height=0.9)
#_, _, peak_lefts, peak_rights = peak_widths(y, peak_centers,rel_height=0.85)
y2 = lininterp_baseline_subtract(y, x, peak_lefts, peak_rights)

# plot to visually inspect baseline subtraction
plt.plot(x, y, '--', label='orig')
plt.plot(x, y2, label='orig - bls')
for i in range(len(peak_centers)):
    plt.plot((x[int(peak_lefts[i] // 1)],
              x[int(peak_rights[i] // 1)]), (y[peak_centers[i]],) * 2, 'b-|')
plt.legend()
plt.show()

y = y2

areas = integrate_peaks(y, peak_lefts, peak_rights)
species = assign_peaks(x, peak_centers, peak_lefts, peak_rights)

wareas = [areas[i] / sensitivity_factors[species[i]]
          for i in range(len(areas))]
total_warea = sum(wareas)

print('Species\tcomposition(mol%)')
for i in range(len(species)):
    comp = 100 * wareas[i] / total_warea
    if comp > 0.1:
        print(f'{species[i]}\t{comp:.2f}')
