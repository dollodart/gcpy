import pathlib as pl
import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.sparse.linalg import spsolve
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


def baseline_als(y, lam, p, niter=10):
    """
    Baseline subtracts by asymmetric least squares.

    Inputs:
      y: vector of data
      lmbd: weighting factor for smoothness
      p: weighting factor for asymmetry
    Outputs:
      z: Least squares solution to Z*w*y
    Locals:
      L: Number of data points in spectrum
      D: Second order approximation of the second derivative
      w: Asymmetric weighting factors
      W: Diagonal matrix of asymmetric weighting factors
      Z: Sum of asymmetric weighting factor and weighted smoothness matrix
    """
    L = len(y)
    D = sparse.diags([1, -2, 1], [0, -1, -2], shape=(L, L - 2))
    w = np.ones(L)
    for i in range(niter):
        W = sparse.spdiags(w, 0, L, L)
        Z = W + lam * D.dot(D.transpose())
        z = spsolve(Z, w * y)
        w = p * (y > z) + (1 - p) * (y < z)
    return z


def find_baseline_inspection(y, l_fineness=8, p_fineness=3, log10l_lims=[
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
            yn = baseline_als(y, lmbd, p)
            e = ((y - yn)**2).mean()
            y2 = np.roll(yn, -2)
            y1 = np.roll(yn, -1)
            s = yn - 2 * y1 + y2
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


def find_baseline_quant(y, l_fineness=8, p_fineness=3,
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
            yn = baseline_als(y, lmbd, p)
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
