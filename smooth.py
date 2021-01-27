import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve, splu


def als_smooth(y, lam, p):
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
    DDl = lam * DD  # weighted smoothness matrix

    w = np.ones(L)  # initialize symmetric least squares
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

    [1] P. H. C. Eilers, "A perfect smoother",
    Anal. Chem. 2003, (75), 3631-3636

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
