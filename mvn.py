#  given covariance matrix and list of indices, tell me what the 
import numpy as _N

def conditionalPDF(mn, cov, mrgidx, cond_on_this=None):
    """
    cov    covariance matrix 
    mrgidx   list of indices to marginalize over  

    #  For given values of 
    """
    dim = cov.shape[0]
    
    allind = _N.arange(dim)
    umrgidx= _N.setdiff1d(allind, mrgidx)

    dimMO= mrgidx.shape[0]
    dimNM = umrgidx.shape[0]
    reord= _N.empty((dim, dim))

    A    = cov[_N.ix_(umrgidx,  umrgidx)]
    B    = cov[_N.ix_(mrgidx,    mrgidx)]
    C    = cov[_N.ix_(mrgidx,   umrgidx)]

    reord[0:dimNM, 0:dimNM] = A
    reord[dimNM:,  dimNM:]  = B
    reord[0:dimNM, dimNM:]      = C.T
    reord[dimNM:,  0:dimNM]     = C

    if cond_on_this is None:
        return A, B, C
    else:
        Bi = _N.linalg.inv(B)
        us = _N.dot(_N.dot(C.T, Bi), cond_on_this)
        return us, A - _N.dot(C.T, _N.dot(Bi, C))

def marginalPDF(mn, cov, mrgidx):
    """
    cov    covariance matrix 
    mrgidx   list of indices to marginalize over  

    #  For given values of 
    """
    dim = cov.shape[0]
    
    allind = _N.arange(dim)
    umrgidx= _N.setdiff1d(allind, mrgidx)

    dimMO= mrgidx.shape[0]
    dimNM = umrgidx.shape[0]
    reord= _N.empty((dim, dim))

    A    = cov[_N.ix_(umrgidx,  umrgidx)]
    B    = cov[_N.ix_(mrgidx,    mrgidx)]
    C    = cov[_N.ix_(mrgidx,   umrgidx)]

    return mn[umrgidx], A

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


