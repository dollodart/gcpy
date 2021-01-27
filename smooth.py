import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve, splu


def als_smooth(y, lam, p):
    """
    Asymmetric least squares smoother.

    Reference: Paul H. C. Eilers, Hans F.M. Boelens.
      "Baseline Correction with Asymmetric Least Squares Smoothing"

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

def whittaker_smooth(y, lmbd):

    """
    Whittaker smoothing algorithm.
    Second order differences used.

    Reference:
    Paul H. C. Eilers, "A perfect smoother",
    Anal. Chem. 2003, (75), 3631-3636


    Inputs:
      y: data vector
      lmbd: smoothing parameter 
    Outputs:
      z: smoothed data
    """

    L = len(y)
    E = sparse.eye(L, format='csc')
    D = sparse.diags([1, -2, 1], [0, -1, -2], shape=(L, L))
    Z = E + lmbd * D.transpose() @ D
    z = spsolve(Z, y)
    return z
