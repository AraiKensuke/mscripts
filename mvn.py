#  given covariance matrix and list of indices, tell me what the 
import numpy as _N


def marginalPDF(cov, lndx, cond_on_this=None):
    """
    cov    covariance matrix 
    lndx   list of indices to marginalize over

    #  For given values of 
    """
    dim = cov.shape[0]
    
    allind = _N.arange(dim)
    ulndx= _N.setdiff1d(allind, lndx)

    dimMO= lndx.shape[0]
    dimM = ulndx.shape[0]
    reord= _N.empty((dim, dim))
    
    A    = cov[_N.ix_(ulndx,  ulndx)]
    B    = cov[_N.ix_(lndx,    lndx)]
    C    = cov[_N.ix_(lndx,   ulndx)]

    reord[0:dimM, 0:dimM] = A
    reord[dimM:,  dimM:]  = B
    B[0:dimM, dimM:]      = C.T
    B[dimM:,  0:dimM]     = C

    if cond_on_this is None:
        return A, B, C
    else:
        Bi = _N.linalg.inv(B)
        us = _N.dot(_N.dot(C.T, Bi), cond_on_this)
        return A, B, C, us, A - _N.dot(C.T, _N.dot(Bi, C))

def practice():    
    #  remove rows and columns from matrix

    A = _N.array([[ 1, -1, -2, -3, -4, -5],
                  [-1,  2,  0,  0,  0,  0],
                  [-2,  0,  3,  0,  0,  0],
                  [-3,  0,  0,  4,  0,  0],
                  [-4,  0,  0,  0,  5,  0],
                  [-5,  0,  0,  0,  0,  6]])

    lndx = _N.array([2, 4])
    ulndx= _N.array([0, 1, 3, 5])

    B = _N.empty((6, 6))
    
    B[0:2, 0:2] = A[_N.ix_(lndx,  lndx)]
    B[2:,  2:]  = A[_N.ix_(ulndx, ulndx)]
    B[0:2, 2:]  = A[_N.ix_(lndx, ulndx)]
    B[2:,  0:2] = A[_N.ix_(ulndx, lndx)]


