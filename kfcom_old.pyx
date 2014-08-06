import numpy as _N
import kfARpnlty as _kfARP

"""
c functions
"""
cdef extern from "math.h":
    double sqrt(double)

#####  common routines for KF and FFBS

#######  Smoothing
def Smooth(_d):
    _d.s_x[_d.N] = _d.f_x[_d.N]
    _d.s_V[_d.N] = _d.f_V[_d.N]

    for n in xrange(_d.N - 1, -1, -1):
        _d.J[n] = _N.dot(_d.f_V[n], _N.dot(_d.F.T, _d.p_Vi[n + 1]))
        _d.s_x[n] = _d.f_x[n] + _N.dot(_d.J[n], _d.s_x[n + 1] - _d.p_x[n + 1])
        #  dgemm has problems handling transpose
        _d.s_V[n] = _d.f_V[n] + _N.dot(_d.J[n], _N.dot(_d.s_V[n + 1] - _d.p_V[n + 1], _d.J[n].T))


##############################################
#################  FFBS    sampling backwarwds

#     return smX
def BS(_d, smXN=None, skip=1):
    GQGT = _N.dot(_d.G, _d.G.T) * _d.Q

    smX = _N.empty((_d.N+1, _d.k))   #  where to store our samples
    if smXN == None:
        smX[_d.N] = _N.random.multivariate_normal(_d.f_x[_d.N,:,0], _d.f_V[_d.N], size=1)
    else:
        smX[_d.N] = smXN[:]

    I        = _N.identity(_d.k)

    for n in xrange(_d.N - 1, -1, -1):
        PFT  =_N.dot(_d.f_V[n], _d.F.T)
        A    =_N.dot(PFT, _N.linalg.inv(_N.dot(_d.F, PFT) + GQGT))
        IAF  = I - _N.dot(A, _d.F)
        PtN  = _N.dot(IAF, _d.f_V[n])
        xm   = _N.dot(IAF, _d.f_x[n, :, 0]) + _N.dot(A, smX[n+1])
        smX[n] = _N.random.multivariate_normal(xm, PtN, size=1)
#        print PtN

    return smX

def BSvec(_d, smXN=None):
    GQGT = _N.dot(_d.G, _d.G.T) * _d.Q

    smX = _N.empty((_d.N+1, _d.k))   #  where to store our samples
    if smXN == None:
        smX[_d.N] = _N.random.multivariate_normal(_d.f_x[_d.N,:,0], _d.f_V[_d.N], size=1)
    else:
        smX[_d.N] = smXN[:]

    fFT    = _N.dot(_d.f_V, _d.F.T)
    #  sum F_{il} V_{lm} F_{mj}
    FfFT   = _N.einsum("il,nlj->nij", _d.F, fFT)
    iv     = _N.linalg.inv(FfFT + _N.tile(GQGT, (_d.N+1,1,1)))
    A      = _N.einsum("nik,nkj->nij", fFT, iv)
    INAF   = _d.IkN - _N.dot(A, _d.F)
    PtN    = _N.einsum("nik,nkj->nij", INAF, _d.f_V)  #  covarainces
    ##  NOW multivariate normal
    mvn1   = _N.random.multivariate_normal(_N.zeros(_d.k), _d.Ik, size=(_d.N+1))
    S,V,D  = _N.linalg.svd(PtN)
    Vs     = _N.sqrt(V)
    VsRn2  =  Vs*mvn1
    zrmn   = _N.einsum("njk,nk->nj", S, VsRn2).reshape((_d.N+1), _d.k, 1)

    #  out of order calculation.  one of the terms can be calculated
    INAFfx = _N.einsum("nj,nj->n", INAF[:, _d.k-1], _d.f_x[:, :, 0])
    last   = _N.zeros(_d.k)
    last[_d.k-1] = 1

    #  temp storage
    Asx = _N.empty(_d.k)
    INAFfxpAsx = _N.empty(_d.k)
    for n in xrange(_d.N - 1, -1, -1):
        #  (I - AF) x
        #smX[n] = zrmn[n, :, 0] + INAFfx[n]*last + _N.dot(A[n], smX[n+1])
        _N.dot(A[n], smX[n+1], out=Asx)
        _N.add(INAFfx[n]*last, Asx, out=INAFfxpAsx)
        _N.add(zrmn[n, :, 0], INAFfxpAsx, out=smX[n])

    return smX

def BS1(_d):
    nrands = _N.random.randn(_d.N+1)
    smX = _N.empty(_d.N+1)
    smX[_d.N] = _d.f_x[_d.N,0,0] + sqrt(_d.f_V[_d.N,0,0]) * nrands[_d.N]
    F0   = _d.F[0,0]
    F02  = F0*F0

    ##########  SMOOTHING step
    GQGT = _d.Q * _d.G[0,0]* _d.G[0,0]
    iGQGT= 1./GQGT

    FTiGQGTF = iGQGT*F02
    for n in xrange(_d.N - 1, -1, -1):
        Sig  = 1/(FTiGQGTF + 1/(_d.f_V[n,0,0]))
        p1   = _d.f_V[n,0,0]* F0
        p2   = 1/(F02*_d.f_V[n,0,0] +GQGT)
        p3   = smX[n+1] - F0*_d.f_x[n,0,0]
        M    = _d.f_x[n,0,0] + p1* p2* p3
        smX[n]= M + sqrt(Sig) * nrands[n]
    return smX
