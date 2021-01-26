import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve, splu

def als_smoother(y, lam, p):
    """
    Asymmetric least squares smoother.

    Inputs:
      y: vector of data
      lmbd: weighting factor for smoothness
      p: weighting factor for asymmetry
         p should always be less than 1, usually less than 0.1
         for positive spectra
         those points greater than the smoothed line are weighted little
         those points lower than the smoothed line are weighted much
    Outputs:
      z: Least squares solution to Z*w*y
         the smoothed spectrum
    Locals:
      L: Number of data points in spectrum
      D: Second order approximation of the second derivative
      w: Asymmetric weighting factors
      W: Diagonal matrix of asymmetric weighting factors
      Z: Linear combination of asymmetric weighting factors
         and weighted smoothness matrix
    """
    L = len(y)

    D = sparse.diags([1, -2, 1], [0, -1, -2], shape=(L, L))
    DD = D.transpose() @ D 
    DDl = lam * DD # weighted smoothness matrix

    w = np.ones(L) # initialize symmetric least squares
    w0 = np.zeros(L) 

    accumulator = 0
    while abs(w - w0).sum() > 0.1 and accumulator < 200:
        w0 = w
        accumulator += 1

        W = sparse.spdiags(w, 0, L, L)
        Z = W + DDl
        z = spsolve(Z, w * y)
        w = p * (y > z) + (1 - p) * (y < z)

    if accumulator >= 200:
        print('warning not converged in 200 iterations')
    return z


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
            yn = als_smoother(y, lmbd, p)
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
            yn = als_smoother(y, lmbd, p)
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

# NOTICE: Code below is copyrighted
# The Whittaker smoother is similar to, but not the same, as the ALS smoother

# Copyright M. H. V. Werts, 2017
#
# martinus point werts à ens-rennes point fr
#
# This software is a computer program whose purpose is to smooth noisy data.
#
# This software is governed by the CeCILL-B license under French law and
# abiding by the rules of distribution of free software.  You can  use,
# modify and/ or redistribute the software under the terms of the CeCILL-B
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and,  more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL-B license and that you accept its terms.
"""
WHITTAKER-EILERS SMOOTHER in Python 3 using numpy and scipy

based on the work by Eilers [1].
    [1] P. H. C. Eilers, "A perfect smoother",
        Anal. Chem. 2003, (75), 3631-3636
coded by M. H. V. Werts (CNRS, France)
tested on Anaconda 64-bit (Python 3.6.4, numpy 1.14.0, scipy 1.0.0)

Read the license text at the end of this file before using this software.

Warm thanks go to Simon Bordeyne who pioneered a first (non-sparse) version
of the smoother in Python.
"""

def speyediff(N, d, format='csc'):
    """
    (utility function)
    Construct a d-th order sparse difference matrix based on
    an initial N x N identity matrix

    Final matrix (N-d) x N
    """

    assert not (d < 0), "d must be non negative"
    shape = (N - d, N)
    diagonals = np.zeros(2 * d + 1)
    diagonals[d] = 1.
    for i in range(d):
        diff = diagonals[:-1] - diagonals[1:]
        diagonals = diff
    offsets = np.arange(d + 1)
    spmat = sparse.diags(diagonals, offsets, shape, format=format)
    return spmat


def whittaker_smooth(y, lmbd, d=2):
    """
    Implementation of the Whittaker smoothing algorithm,
    based on the work by Eilers [1].

    [1] P. H. C. Eilers, "A perfect smoother", Anal. Chem. 2003, (75), 3631-3636

    The larger 'lmbd', the smoother the data.
    For smoothing of a complete data series, sampled at equal intervals

    This implementation uses sparse matrices enabling high-speed processing
    of large input vectors

    ---------

    Arguments :

    y       : vector containing raw data
    lmbd    : parameter for the smoothing algorithm (roughness penalty)
    d       : order of the smoothing

    ---------

    Returns :

    z       : vector of the smoothed data.
    """

    m = len(y)
    E = sparse.eye(m, format='csc')
    D = speyediff(m, d, format='csc')
    coefmat = E + lmbd * D.conj().T.dot(D)
    z = splu(coefmat).solve(y)
    return z


