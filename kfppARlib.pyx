import numpy as _N
cimport numpy as _N   #  need to add include path  /Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/numpy/core/include

import scipy.optimize as _sco
import kfpp as _kfpp
import kfcom as _kfcom

DTYP   = _N.float64
ctypedef _N.float64_t DTYP_t

"""
c functions
"""
cdef extern from "math.h":
    double exp(double)
    double sqrt(double)
    double log(double)
    double abs(double)

#cdef double dt = 0.0005
cdef double dt = 0.001


########################   KF
def armdl_KF_1itr(_d, approx=0):   #  approximation
    if approx == 0:          FF(_d)
    elif approx == 1:        FFa(_d)
    elif approx == 2:        FFaa(_d)
    _kfcom.Smooth(_d)

########################   FFBS
def armdl_FFBS_1itr(_d, approx=0, samples=1):   #  approximation
    if approx == 0:          FF(_d)
    elif approx == 1:        FFa(_d)
    elif approx == 2:        FF(_d)
    if samples > 1:
        smpls = []
        for s in xrange(samples):
            smpls.append(_kfcom.BS(_d))
    else:
        if _d.k == 1:
            smpls = _kfcom.BS1(_d)
        else:
            smpls = _kfcom.BS(_d)
    return smpls

##############################
##############################
def FF(_d):
    GQGT= _d.Q * _N.dot(_d.G, _d.G.T)
    bbT = _N.dot(_d.beta, _d.beta.T)
    cdef int n
    zrs = _N.empty(_d.N)
    B0  = _d.beta[0, 0]
    for n from 1 <= n < _d.N + 1:
        _d.p_x[n] = _N.dot(_d.F, _d.f_x[n - 1])
        _d.p_V[n] = _N.dot(_d.F, _N.dot(_d.f_V[n - 1], _d.F.T)) + GQGT
        #  set _d.f_x
        if _d.k == 1:
            _d.p_Vi[n, 0, 0] = 1 / _d.p_V[n, 0, 0]   #  invert
        else:
            _d.p_Vi[n] = _N.linalg.inv(_d.p_V[n])   #  invert

        sol    = _sco.root(_kfpp.f_and_df, _d.p_x[n, :, 0], jac=True, args=(_d.k, _d.p_x[n], _d.p_V[n], _d.p_Vi[n], _d.beta, _d.u, _d.dN[n], bbT), tol=1e-5)
        _d.f_x[n, :, 0] = sol.x

        ef = exp(_d.u + B0*_d.f_x[n,0,0])*dt
        trm1 = ef / ((1+ef)*(1+ef))
        if _d.k == 1:
            _d.f_V[n, 0, 0] = 1 / (_d.p_Vi[n, 0, 0] + bbT * trm1)
        else:
            _d.f_V[n] = _N.linalg.inv(_d.p_Vi[n] + bbT * trm1)


#  replace Laplace method only
##############################
##############################
def FFa(_d):
    GQGT= _d.Q * _N.dot(_d.G, _d.G.T)
    bbT = _N.dot(_d.beta, _d.beta.T)
    cdef int n
    zrs = _N.empty(_d.N)
    B0  = _d.beta[0, 0]

    lstSpk = _d.N

    for n from 1 <= n < _d.N + 1:
        _d.p_x[n] = _N.dot(_d.F, _d.f_x[n - 1])
        _d.p_V[n] = _N.dot(_d.F, _N.dot(_d.f_V[n - 1], _d.F.T)) + GQGT
        #  set _d.f_x

        if _d.k == 1:
            _d.p_Vi[n, 0, 0] = 1 / _d.p_V[n, 0, 0]   #  invert
        else:
            _d.p_Vi[n] = _N.linalg.inv(_d.p_V[n])   #  invert

        ######  APPROX LAPLACE #######
        if (n > 40) and ((_d.dN[n] == 0) and (n - lstSpk > _d.tau)):  #APPROX
            if _d.k == 1:
                _d.f_V[n, 0, 0] = 1 / (_d.p_Vi[n, 0, 0] + bbT * exp(_d.u + B0*_d.p_x[n, 0, 0])*dt)
            else:
                _d.f_V[n] = _N.linalg.inv(_d.p_Vi[n] + bbT * exp(_d.u + B0*_d.p_x[n, 0, 0])*dt)
            _d.f_x[n] = _d.p_x[n] + _N.dot(_d.f_V[n], _d.beta) * (_d.dN[n] - exp(_d.u + B0*_d.p_x[n, 0, 0])*dt)
        ######  EXACT LAPLACE #######
        else:
            if _d.dN[n] == 1:
                lstSpk = n

            sol    = _sco.root(_kfpp.f_and_df, _d.p_x[n, :, 0], jac=True, args=(_d.k, _d.p_x[n], _d.p_V[n], _d.p_Vi[n], _d.beta, _d.u, _d.dN[n], bbT), tol=1e-5)
            _d.f_x[n, :, 0] = sol.x
            if _d.k == 1:
                _d.f_V[n, 0, 0] = 1/(_d.p_Vi[n, 0, 0] + B0*B0 * exp(_d.u + B0 * _d.f_x[n, 0, 0])*dt)
            else:
                _d.f_V[n] = _N.linalg.inv(_d.p_Vi[n] + bbT * exp(_d.u + B0 * _d.f_x[n, 0, 0])*dt)

#  For k == 1, this method not necessary, because inverse costs nothing
#  For some reason, when KF_1itr is called before KFaa_1itr, aa works fine
#  when aa called first, after first
########################
########################   Forward filter
def FFaa(_d):   #  approximation
    print "inv. approximated both"
    GQGT= _d.Q * _N.dot(_d.G, _d.G.T)
    bbT = _N.dot(_d.beta, _d.beta.T)
    cdef int n
    zrs = _N.empty(_d.N)

    lstSpk = _d.N

    B0     = _d.beta[0, 0]

    Nx  = 50
    ##########  First 50 steps, exact Laplace
    for n from 1 <= n < Nx:
        _d.p_x[n] = _N.dot(_d.F, _d.f_x[n - 1])
        _d.p_V[n] = _N.dot(_d.F, _N.dot(_d.f_V[n - 1], _d.F.T)) + GQGT
        #  set _d.f_x
        if _d.k == 1:
            _d.p_Vi[n, 0, 0] = 1 / _d.p_V[n, 0, 0]   #  invert
        else:
            _d.p_Vi[n] = _N.linalg.inv(_d.p_V[n])   #  invert

        sol    = _sco.root(_kfpp.f_and_df, _d.p_x[n, :, 0], jac=True, args=(_d.k, _d.p_x[n], _d.p_V[n], _d.p_Vi[n], _d.beta, _d.u, _d.dN[n], bbT), tol=1e-5)
        _d.f_x[n, :, 0] = sol.x
        ef = exp(_d.u + B0 * _d.f_x[n,0,0])*dt
        if _d.k == 1:
            _d.f_V[n, 0, 0] = 1/(_d.p_Vi[n, 0, 0] + B0*B0 * ef / ((1+ef)*(1+ef)))
        else:
            _d.f_V[n] = _N.linalg.inv(_d.p_Vi[n] + bbT * ef / ((1+ef)*(1+ef)))

    if _d.k > 1:
        ep   = exp(_d.u + B0*_d.p_x[Nx, 0, 0])*dt
        Vinmnm = _d.p_Vi[Nx-1] + bbT * ep / ((1+ep)*(1+ep))
    ##########  Decide exact or approximate Laplace
    for n from Nx <= n < _d.N + 1:
        _d.p_x[n] = _N.dot(_d.F, _d.f_x[n - 1])
        _d.p_V[n] = _N.dot(_d.F, _N.dot(_d.f_V[n - 1], _d.F.T)) + GQGT
        ######  APPROX LAPLACE #######
        if (_d.dN[n] == 0) and (n - lstSpk > _d.tau):
            ep   = exp(_d.u + B0*_d.p_x[n, 0, 0])*dt
            if _d.k == 1:  #  Simple scalar inversion
                _d.p_Vi[n, 0, 0] = 1/_d.p_V[n, 0, 0]   #  invert
                _d.f_V[n] = 1 / (_d.p_Vi[n, 0, 0] + bbT * ep/((1+ep)*(1+ep)))
            else:          #  Matrix inversion approx
                #  FIRST approximate inv mat:   find _d.p_Vi[n]
                e   = _N.dot(_d.p_Vi[n - 1], _d.p_V[n] - _d.p_V[n - 1])
                _d.p_Vi[n] = _N.dot(_d.Ik - e, _d.p_Vi[n - 1])
                #  SECOND approximate inv mat:  find _d.f_V[n]
                Vinn = _d.p_Vi[n] + bbT * ep / ((1+ep)*(1+ep))   # inverse filter covariance
                e   = _N.dot(_d.f_V[n - 1], Vinn - Vinmnm)
                _d.f_V[n] = _N.dot(_d.Ik - e, _d.f_V[n - 1])
            #  we first need _d.f_V[n]

            #  replacement for Laplace
            _d.f_x[n] = _d.p_x[n] + _N.dot(_d.f_V[n], _d.beta) * (_d.dN[n] - ep / (1 + ep))
        ######  EXACT LAPLACE #######
        else:
            if _d.dN[n] == 1:  lstSpk = n

            if _d.k == 1:    _d.p_Vi[n, 0, 0] = 1/_d.p_V[n, 0, 0]
            else:            _d.p_Vi[n] = _N.linalg.inv(_d.p_V[n])

            sol    = _sco.root(_kfpp.f_and_df, _d.p_x[n, :, 0], jac=True, args=(_d.k, _d.p_x[n], _d.p_V[n], _d.p_Vi[n], _d.beta, _d.u, _d.dN[n], bbT), tol=1e-5)
            _d.f_x[n, :, 0] = sol.x
            ef   = exp(_d.u + B0*_d.f_x[n, 0, 0])*dt
            if _d.k == 1:
                _d.f_V[n, 0, 0] = 1/(_d.p_Vi[n, 0, 0] + bbT*ef/((1+ef)*(1+ef)))
            else:
                Vinn = _d.p_Vi[n] + bbT * ef / ((1+ef)*(1+ef))   #  evaluate here, used for approximation in next step
                _d.f_V[n] = _N.linalg.inv(Vinn)

        if _d.k != 1:
            Vinmnm = Vinn

def armdl_logL(_d):
    cdef unsigned int n
    logL1 = 0
    logL2 = 0
    for n from 0 <= n < _d.N + 1:
        pe= exp(_d.u + _d.beta[0, 0] * _d.s_x[n, 0, 0])*dt
        if _d.dN[n] == 1:
            logL1 += _N.log(pe / (1 + pe))
        else:
            logL1 += _N.log(1  / (1 + pe))

    return logL1

