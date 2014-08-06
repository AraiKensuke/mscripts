import numpy as _N
cimport numpy as _N   #  need to add include path  /Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/numpy/core/include
import kfcom as _kfcom

import scipy.optimize as _sco

import kfARpnlty as _kfARP
import time as _tm

import warnings
warnings.filterwarnings("error")

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
    FF(_d)
    _kfcom.Smooth(_d)

########################   FFBS
def armdl_FFBS_1itr(_d, m=0, offset=0, samples=1, ffast=False, fast=False, smXN=None):   #  approximation
   ##########  FF
    if _d.k == 1:
        FF1dv(_d, offset=offset)
    else:
        FFdv(_d, m, fast=fast)
   ##########  BS
    if _d.k == 1:
        smpls = _kfcom.BS1(_d)
    else:
        #try:
        smXN = _N.random.multivariate_normal(_d.f_x[m, _d.N,:,0], _d.f_V[m, _d.N], size=1)

        if ffast:
            smpls = _kfcom.BSvec(_d, m=m, smXN=smXN)
        else:
            smpls = _kfcom.BS(_d, m=m, smXN=smXN)

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

def FFdv(_d, m=0, offset=0, fast=False):   #  approximate KF    #  k==1,dynamic variance
    GQGT    = _N.dot(_d.G, _d.G.T)*_d.Q[m]
    #  do this until p_V has settled into stable values
    
    k     = _d.k
    px    = _d.p_x[m]
    pV    = _d.p_V[m]
    fx    = _d.f_x[m]
    fV    = _d.f_V[m]
    Rv    = _d.Rv[m]
    K     = _d.K[m]

    """
    temporary storage
    """
    Hpx   = _N.empty((1, 1))
    #KyoHpx= _N.empty((k, 1))
    KH    = _N.empty((k, k))    
    IKH   = _N.empty((k, k))
    VFT   = _N.empty((k, k))
    FVFT  = _N.empty((k, k))
    yo    = _d.y[m]

    if not fast:
        for n from 1 <= n < _d.N + 1:
            px[n] = _N.dot(_d.F, fx[n - 1])
            pV[n] = _N.dot(_d.F, _N.dot(fV[n - 1], _d.F.T)) + GQGT
            mat  = 1 / (pV[n, 0, 0] + Rv[n])
            K[n] = _N.dot(pV[n], mat*_d.H.T)   #  vector
            fx[n]    = px[n] + _N.dot(K[n], (yo[n] - _N.dot(_d.H, px[n])))
            fV[n] = _N.dot(_d.Ik - _N.dot(K[n], _d.H), pV[n])
    else:
        for n from 1 <= n < _d.N + 1:
            _N.dot(_d.F, fx[n - 1], out=px[n])
            _N.dot(fV[n - 1], _d.F.T, out=VFT)
            _N.dot(_d.F, VFT, out=FVFT)
            _N.add(FVFT, GQGT, out=pV[n])
            mat  = 1 / (pV[n, 0, 0] + Rv[n])
            _N.dot(pV[n], mat*_d.H.T, out=K[n])   #  vector

            # px + K(y - o - Hpx)  K column vec, (y-o-Hpx) is scalar
            _N.dot(_d.H, px[n], out=Hpx)
            KyoHpx = K[n]* (yo[n] - Hpx[0, 0])
            _N.add(px[n], KyoHpx, out=fx[n])
                                  
            # (I - KH)pV
            _N.dot(K[n], _d.H, out=KH)
            _N.subtract(_d.Ik, KH, out=IKH)
            _N.dot(IKH, pV[n], out=fV[n])

def FF1dv(_d, offset=0):   #  approximate KF    #  k==1,dynamic variance
    GQGT    = _d.G[0,0]*_d.G[0, 0] * _d.Q
    k     = _d.k
    px    = _d.p_x
    pV    = _d.p_V
    fx    = _d.f_x
    fV    = _d.f_V
    Rv    = _d.Rv
    K     = _d.K

    #  do this until p_V has settled into stable values

    for n from 1 <= n < _d.N + 1:
        px[n,0,0] = _d.F[0,0] * fx[n - 1,0,0]
#        pV[n,0,0] = _d.F[0,0] * fV[n - 1,0,0] * _d.F.T[0,0] + GQGT
        pV[n,0,0] = _d.F[0,0] * fV[n - 1,0,0] * _d.F[0,0] + GQGT
        #_d.p_Vi[n,0,0] = 1/pV[n,0,0]

#        mat  = 1 / (_d.H[0,0]*pV[n,0,0]*_d.H[0,0] + Rv[n])
        mat  = 1 / (pV[n,0,0] + Rv[n])
#        K[n,0,0] = pV[n]*_d.H[0,0]*mat
        K[n,0,0] = pV[n,0,0]*mat
#        fx[n,0,0]    = px[n,0,0] + K[n,0,0]*(_d.y[n] - offset[n] - _d.H[0,0]* px[n,0,0])
#        fx[n,0,0]    = px[n,0,0] + K[n,0,0]*(_d.y[n] - _d.H[0,0]* px[n,0,0])
        fx[n,0,0]    = px[n,0,0] + K[n,0,0]*(_d.y[n] - px[n,0,0])
#        fV[n,0,0] = (1 - K[n,0,0]* _d.H[0,0])* pV[n,0,0]
        fV[n,0,0] = (1 - K[n,0,0])* pV[n,0,0]


