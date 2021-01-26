from gcpy.decode import read_data
from gcpy.smooth import als_smooth
import matplotlib.pyplot as plt
import numpy as np

def als_smoothing_inspection(y, l_fineness=8, p_fineness=3, log10l_lims=[
                             2, 9], log10p_lims=[-3, 0]):
    """

    Provides a set of plots for inspection of optimal parameters.
    Default values make the search space suggested in the paper.

    """

    lmbd_range = 10.**np.linspace(log10l_lims[0], log10l_lims[1], l_fineness)
    p_range = 10.**np.linspace(log10p_lims[0], log10p_lims[1], p_fineness)
    fig, axs = plt.subplots(nrows=len(lmbd_range), ncols=len(p_range))
    for c, lmbd in enumerate(lmbd_range):
        for d, p in enumerate(p_range):
            yn = als_smooth(y, lmbd, p)
            e = ((y - yn)**2).mean()
            #y2 = np.roll(yn, -2)
            #y1 = np.roll(yn, -1)
            #s = yn - 2 * y1 + y2
            axs[c, d].plot(y)
            axs[c, d].plot(yn)
            axs[c, d].set_xticks([])
            axs[c, d].set_yticks([])
    axs[0, 0].set_title(
        r'$\log_{{10}} p=[{0},{1}],{2},\rightarrow$'.format(*log10p_lims, p_fineness))
    axs[0, -
        1].set_title(r'$\log_{{10}} \lambda=[{0},{1}],{2},\downarrow$'.format(*
                                                                              log10l_lims, l_fineness))
    # updownarrow does not have unicode equiv
    return fig, axs


def als_smoothing_quantification(y, l_fineness=8, p_fineness=3,
                        log10l_lims=[2, 9], log10p_lims=[-3, 0]):
    """

    Makes plots of the root mean square error and the smoothness (as
    measured to second order approximation).  Assumes unit spacing in x
    variable.

    """

    mse = []
    smo = []
    handles = []
    lmbd_range = 10.**np.linspace(log1 - l_lims[0], log10l_lims[1], l_fineness)
    p_range = 10.**np.linspace(log10p_lims[0], log10p_lims[1], p_fineness)

    for c, lmbd in enumerate(lmbd_range):
        for d, p in enumerate(p_range):
            yn = als_smooth(y, lmbd, p)
            e = ((y - yn)**2).mean()
            y2 = np.roll(yn, -2)
            y1 = np.roll(yn, -1)
            s = yn - 2 * y1 + y2
            mse.append(e)
            s = ((s**2 * yn)**2).mean()
            smo.append(s)

    mse = np.array(mse).reshape(len(lmbd_range), len(p_range))
    smo = np.array(smo).reshape(len(lmbd_range), len(p_range))
    lmbd, p = np.meshgrid(lmbd_range, p_range, indexing='ij')

    fig, axs = plt.subplots(nrows=1, ncols=2)
    CS = axs[0].contour(np.log10(lmbd), np.log10(p), mse)
    axs[0].clabel(CS, inline=True, fmt='{0:.2E}'.format)
    axs[0].set_xlabel(r'$\log_{10}\lambda$')
    axs[0].set_ylabel(r'$\log_{10}p$')
    axs[0].set_title(r'$\langle (y-y_s)^2\rangle$')

    CS = axs[1].contour(np.log10(lmbd), np.log10(p), smo)
    axs[1].clabel(CS, inline=True, fmt='{0:.2E}'.format)
    axs[1].set_xlabel(r'$\log_{10}\lambda$')
    axs[1].set_ylabel(r'$\log_{10}p$')
    axs[1].set_title(
        r'$\langle \langle \Delta^2 z)^2 \rangle$ where $\Delta\equiv \frac{\delta^2 y}{\delta x^2}$')

    return fig, axs

# read in data
data = read_data('../data/baseline-and-noisy-test-data.CH')
x = data['time']
y = data['tic']
als_smoothing_inspection(y, l_fineness=5, p_fineness=3, log10l_lims=[1, 5], log10p_lims=[-1, -0.301])
plt.show()
