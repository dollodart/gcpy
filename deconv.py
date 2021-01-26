import matplotlib.pyplot as plt
import numpy as np
import numpy.random as r
from lmfit import minimize, Parameters
from rampy_funcs import *
from encode import encode
from decode import read_data
from scipy.integrate import trapz
from scipy.signal import find_peaks
from bls import lininterp_baseline_subtract

x = np.linspace(0, 1, 1000)
ampg = 1
mg = 0.5
sg = 0.05
ampp = 0.8
mp = 0.6
sp = 0.2
lg = 0.2

y1 = gaussian(x, ampg, mg, sg)
y2 = pseudovoigt(x, ampp, mp, sp, lg)
encode('data/deconv.CH', y1 + y2 + r.random(len(y1)) / 10., xmin=0, xmax=1)

data = read_data('data/deconv.CH')
x = data['time']
y = data['tic']
peak_pos, peak_prop = find_peaks(y, width=40, height=1, rel_height=0.99)
y = lininterp_baseline_subtract(
    y, x, peak_prop['left_ips'], peak_prop['right_ips'])

#plt.plot(x, y)
# plt.show()
#import sys; sys.exit()

# least square fit using lmfit


def c1(params, x):
    return gaussian(x,
                    params['ampg'],
                    params['mg'],
                    params['sg'])


def c2(params, x):
    return pseudovoigt(
        x,
        params['ampp'],
        params['mp'],
        params['sp'],
        params['LG'])


def model(params, x):
    return c1(params, x) + c2(params, x)


def residual(params, x, data):
    return data - model(params, x)


params = Parameters()
params.add('ampg', value=ampg + r.random_sample() * ampg / 4)
params.add('mg', value=mg + r.random_sample() * mg / 4)
params.add('sg', value=sg + r.random_sample() * sg / 4)
params.add('ampp', value=ampp + r.random_sample() * ampp / 4)
params.add('mp', value=mp + r.random_sample() * mp / 4)
params.add('sp', value=sp + r.random_sample() * sp / 4)
params.add('LG', value=lg + r.random_sample() * lg / 4)

out = minimize(residual, params, args=(x, y))

# visualize fit and residuals
fig = plt.figure()
ax1 = fig.add_subplot(2, 1, 1)
ax2 = fig.add_subplot(2, 1, 2)  # ,sharex=True)
ax1.plot(x, model(out.params, x))
ax1.plot(x, c1(out.params, x))
ax1.plot(x, c2(out.params, x))
ax1.plot(x, y)
ax2.plot(x, residual(out.params, x, y))
plt.show()

# quantify independent peak areas
# there are analytical formulas for the integration area of some profiles
c1_area = trapz(c1(out.params, x))
c2_area = trapz(c2(out.params, x))
m_area = trapz(y)
print(
    "C1={0:.1f},C2={1:.1f} mol %".format(
        100 *
        c1_area /
        m_area,
        100 *
        c2_area /
        m_area))
