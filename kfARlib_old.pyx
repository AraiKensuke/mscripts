import numpy as _N
cimport numpy as _N   #  need to add include path  /Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/numpy/core/include
import kfcom as _kfcom

import scipy.optimize as _sco

import kfARpnlty as _kfARP
import time as _tm

"""
c functions
"""
cdef extern from "math.h":
    double exp(double)
    double sqrt(double)
    double log(double)
    double abs(double)
"""
p        AR order
Ftrgt    Ftrgt[0]  noise amp.  Ftrgt[1:]  AR(p) coeffs
f        freq.
f0       bandpass
f1
zr       amp. at band stop
"""

########################   KF
def armdl_KF_1itr(_d, approx=0):   #  approximation
    if approx == 0:          FF(_d)
    elif approx == 1:        FFa(_d)
    _kfcom.Smooth(_d)

########################   FFBS
def armdl_FFBS_1itr(_d, offset=0, samples=1, ffast=False, fast=False, smXN=None):   #  approximation
    if _d.k == 1:
        FF1dv(_d, offset=offset)
    else:
        #t1 = _tm.time()
        FFdv(_d, fast=ffast)
        #t2 = _tm.time()
    if samples > 1:
        smpls = []
        for s in xrange(samples):
            smpls.append(_kfcom.BS(_d))
    else:
        if _d.k == 1:
            smpls = _kfcom.BS1(_d)
        else:
            smXN = _N.random.multivariate_normal(_d.f_x[_d.N,:,0], _d.f_V[_d.N], size=1)
            if fast:
                smpls = _kfcom.BSvec(_d, smXN=smXN)
            else:
                smpls = _kfcom.BS(_d, smXN=smXN)

    return smpls

def FF(_d):   #  approximate KF
    GQGT    = _N.dot(_d.G, _N.dot(_d.Q, _d.G.T))

    cdef int n
    #  do this until p_V has settled into stable values
    for n from 1 <= n < _d.N + 1:
        _d.p_x[n] = _N.dot(_d.F, _d.f_x[n - 1])
        _d.p_V[n] = _N.dot(_d.F, _N.dot(_d.f_V[n - 1], _d.F.T)) + GQGT
        _d.p_Vi[n] = _N.linalg.inv(_d.p_V[n])

        mat  = 1 / (_N.dot(_d.H, _N.dot(_d.p_V[n], _d.H.T)) + _d.R)
        _d.K[n] = _N.dot(_d.p_V[n], _N.dot(_d.H.T, mat))
        _d.f_x[n]    = _d.p_x[n] + _N.dot(_d.K[n], (_d.y[n] - _N.dot(_d.H, _d.p_x[n])))
        _d.f_V[n] = _N.dot(_d.Ik - _N.dot(_d.K[n], _d.H), _d.p_V[n])

def FFdv(_d, offset=0, fast=False):   #  approximate KF    #  k==1,dynamic variance
    GQGT    = _N.dot(_d.G, _N.dot(_d.Q, _d.G.T))
    #  do this until p_V has settled into stable values

    
    """
    temporary storage
    """
    Hpx   = _N.empty((1, 1))
    #KyoHpx= _N.empty((_d.k, 1))
    KH    = _N.empty((_d.k, _d.k))    
    IKH   = _N.empty((_d.k, _d.k))
    VFT   = _N.empty((_d.k, _d.k))
    FVFT  = _N.empty((_d.k, _d.k))
    yo = _d.y - offset
    if not fast:
        for n from 1 <= n < _d.N + 1:
            _d.p_x[n] = _N.dot(_d.F, _d.f_x[n - 1])
            _d.p_V[n] = _N.dot(_d.F, _N.dot(_d.f_V[n - 1], _d.F.T)) + GQGT
            mat  = 1 / (_d.p_V[n, 0, 0] + _d.Rv[n])
            _d.K[n] = _N.dot(_d.p_V[n], mat*_d.H.T)   #  vector
            _d.f_x[n]    = _d.p_x[n] + _N.dot(_d.K[n], (_d.y[n] - offset - _N.dot(_d.H, _d.p_x[n])))
            _d.f_V[n] = _N.dot(_d.Ik - _N.dot(_d.K[n], _d.H), _d.p_V[n])
    else:
        for n from 1 <= n < _d.N + 1:
            _N.dot(_d.F, _d.f_x[n - 1], out=_d.p_x[n])
            _N.dot(_d.f_V[n - 1], _d.F.T, out=VFT)
            _N.dot(_d.F, VFT, out=FVFT)
            _N.add(FVFT, GQGT, out=_d.p_V[n])
            mat  = 1 / (_d.p_V[n, 0, 0] + _d.Rv[n])
            _N.dot(_d.p_V[n], mat*_d.H.T, out=_d.K[n])   #  vector

            # px + K(y - o - Hpx)  K column vec, (y-o-Hpx) is scalar
            _N.dot(_d.H, _d.p_x[n], out=Hpx)
            KyoHpx = _d.K[n]* (yo[n] - Hpx[0, 0])
            _N.add(_d.p_x[n], KyoHpx, out=_d.f_x[n])
                                  
            # (I - KH)pV
            _N.dot(_d.K[n], _d.H, out=KH)
            _N.subtract(_d.Ik, KH, out=IKH)
            _N.dot(IKH, _d.p_V[n], out=_d.f_V[n])

def FF1dv(_d, offset=0):   #  approximate KF    #  k==1,dynamic variance
    GQGT    = _d.G[0,0]*_d.G[0, 0] * _d.Q

    #  do this until p_V has settled into stable values

    for n from 1 <= n < _d.N + 1:
        _d.p_x[n,0,0] = _d.F[0,0] * _d.f_x[n - 1,0,0]
#        _d.p_V[n,0,0] = _d.F[0,0] * _d.f_V[n - 1,0,0] * _d.F.T[0,0] + GQGT
        _d.p_V[n,0,0] = _d.F[0,0] * _d.f_V[n - 1,0,0] * _d.F[0,0] + GQGT
        _d.p_Vi[n,0,0] = 1/_d.p_V[n,0,0]

#        mat  = 1 / (_d.H[0,0]*_d.p_V[n,0,0]*_d.H[0,0] + _d.Rv[n])
        mat  = 1 / (_d.p_V[n,0,0] + _d.Rv[n])
#        _d.K[n,0,0] = _d.p_V[n]*_d.H[0,0]*mat
        _d.K[n,0,0] = _d.p_V[n,0,0]*mat
#        _d.f_x[n,0,0]    = _d.p_x[n,0,0] + _d.K[n,0,0]*(_d.y[n] - offset[n] - _d.H[0,0]* _d.p_x[n,0,0])
#        _d.f_x[n,0,0]    = _d.p_x[n,0,0] + _d.K[n,0,0]*(_d.y[n] - _d.H[0,0]* _d.p_x[n,0,0])
        _d.f_x[n,0,0]    = _d.p_x[n,0,0] + _d.K[n,0,0]*(_d.y[n] - _d.p_x[n,0,0])
#        _d.f_V[n,0,0] = (1 - _d.K[n,0,0]* _d.H[0,0])* _d.p_V[n,0,0]
        _d.f_V[n,0,0] = (1 - _d.K[n,0,0])* _d.p_V[n,0,0]


########################
########################   KF
#cimport cython
#@cython.boundscheck(False)
def FFa(_d):   #  approximate KF
    GQGT    = _N.dot(_d.G, _N.dot(_d.Q, _d.G.T))

    cdef int n
    #  do this until p_V has settled into stable values
    for n from 1 <= n < 40:
        _d.p_x[n] = _N.dot(_d.F, _d.f_x[n - 1])
        _d.p_V[n] = _N.dot(_d.F, _N.dot(_d.f_V[n - 1], _d.F.T)) + GQGT
        _d.p_Vi[n] = _N.linalg.inv(_d.p_V[n])

        mat  = 1 / (_N.dot(_d.H, _N.dot(_d.p_V[n], _d.H.T)) + _d.R)
        _d.K[n] = _N.dot(_d.p_V[n], _N.dot(_d.H.T, mat))
        _d.f_x[n]    = _d.p_x[n] + _N.dot(_d.K[n], (_d.y[n] - _N.dot(_d.H, _d.p_x[n])))
        _d.f_V[n] = _N.dot(_d.Ik - _N.dot(_d.K[n], _d.H), _d.p_V[n])

    ########
    for n from 40 <= n < _d.N + 1:
        _d.p_x[n] = _N.dot(_d.F, _d.f_x[n - 1])
        _d.p_V[n] = _N.dot(_d.F, _N.dot(_d.f_V[n - 1], _d.F.T)) + GQGT
        #  I know _d.p_Vi[n - 1]
        #  e = (Vn|n-1 - Vn-1|n-2) Vin-1|n-2
        #  Vin|n-1 = (1 - e)Vin-1|n-2   #  don't compute matrix inverse
        e   = _N.dot(_d.p_Vi[n - 1], _d.p_V[n] - _d.p_V[n - 1])
        _d.p_Vi[n] = _N.dot(_d.Ik - e, _d.p_Vi[n - 1])
        mat  = 1 / (_N.dot(_d.H, _N.dot(_d.p_V[n], _d.H.T)) + _d.R)
        _d.K[n] = _N.dot(_d.p_V[n], _N.dot(_d.H.T, mat))
        _d.f_x[n]    = _d.p_x[n] + _N.dot(_d.K[n], (_d.y[n] - _N.dot(_d.H, _d.p_x[n])))
        _d.f_V[n] = _N.dot(_d.Ik - _N.dot(_d.K[n], _d.H), _d.p_V[n])

def armdl_logL(_d):
    cdef unsigned int n
#    cdef double logL = -((_d.N + 1)/2.) * log(_N.pi * _d.R)
    cdef double logL = -(_d.N + 1) * log(_N.pi * _d.R)

    for n from 0 <= n < _d.N + 1:
        x = _d.y[n] - _d.H[0, 0] * _d.s_x[n, 0, 0]
        logL -= (x*x) / (2*_d.R)

    return logL

def armdl_logL2(_d):
    cdef unsigned int n
#    cdef double logL = -((_d.N + 1)/2.) * log(_N.pi * _d.R)
    cdef double logL = 0

    H0   = _d.H[0, 0]
    for n from 1 <= n < _d.N + 1:
        logL += log(H0*H0*_d.p_V[n, 0, 0] + _d.R)
        logL += (_d.y[n] - H0*_d.p_x[n, 0, 0])*(H0*H0*_d.p_V[n, 0, 0] + _d.R)*(_d.y[n] - H0*_d.p_x[n, 0, 0])
    return -0.5 * logL
