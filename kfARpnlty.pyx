import numpy as _N
#cimport numpy as _N   #  need to add include path  /Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/numpy/core/include

import scipy.optimize as _sco
import scipy.integrate as _sciint

#DTYP   = _N.float64
#ctypedef _N.float64_t DTYP_t

"""
c functions
"""
cdef extern from "math.h":
    double cos(double)
    double log(double)

def _b(f, F0, FF0, Indx, IIndx, p):
    b = 1
    b -= 2*_N.dot(F0, _N.cos(f*Indx))
    b += _N.sum(FF0*_N.cos(f*IIndx))  #  FF0, IIndx are 2-D
    return b

def _db(f, F0, Indx, p, k):
    db = -2 * cos(6.283185307179586*k*f)
    db += 2*_N.dot(F0, _N.cos(f*(6.283185307179586*k - Indx))) #  k -> 2*_N.pi*k
    return db

def _ddb(f, F0, p, k, l):
    return 2*cos(6.283185307179586*(k-l)*f)

def mdQlJ(Fnz, _d, _e):
    """
    Fnz    low and high AR + noizes.  Need row vector
    """
    #  components 1 + p_L + 1 + p_H  (incl. noise term)
    pq_L = 1 + _d.kL;         pq_H = 1 + _d.kH
    Fnz= Fnz.reshape((pq_L + pq_H, 1))
    #  Fnz[0]    noize   Fnz[1]        ... Fnz[pq_L - 1]       dynamics
    #  Fnz[pq_L] noize   Fnz[pq_L + 1] ... Fnz[pq_L + pqH - 1] dynamics
    q2_L  = Fnz[0, 0];      q2_H  = Fnz[pq_L, 0]
    Fl  = Fnz[1:pq_L, 0];   Fh  = Fnz[pq_L+1:, 0]
    F2l = _N.dot(Fnz[1:pq_L, :], Fnz[1:pq_L, :].T)
    F2h = _N.dot(Fnz[pq_L + 1, :], Fnz[pq_L + 1, :].T)

    q4_L = q2_L*q2_L;    q6_L = q4_L*q2_L
    q4_H = q2_H*q2_H;    q6_H = q4_H*q2_H

    dQmJ = _N.empty(pq_L + pq_H)
    dQmJ[1:pq_L]   = (_e.BL - _N.dot(_e.AL, Fl)) / q2_L    # incomplete
    dQmJ[0]         = -_d.N / (2*q2_L) + (_e.CL - 2*_N.dot(Fl.T, _e.BL) + _N.dot(Fl.T, _N.dot(_e.AL, Fl))) / (2*q4_L) - _e.lmbda*(_e.fL2 - _e.fL1)/q2_L
    #######################
    dQmJ[1+pq_L:]    = (_e.BH - _N.dot(_e.AH, Fh)) / q2_H # incomplete
    dQmJ[pq_L]   = -_d.N / (2*q2_H) + (_e.CH - 2*_N.dot(Fh.T, _e.BH) + _N.dot(Fh.T, _N.dot(_e.AH, Fh))) / (2*q4_H) - _e.lmbda*(_e.fH2 - _e.fH1)/q2_H

    #########  LOW
    for k from 1 <= k < pq_L:
        intgrl = _sciint.quad(intg1, _e.fL1, _e.fL2, args=(Fl, F2l, _e.IxL, _e.IIxL, _d.kL, k), epsrel=_e.epsrel, epsabs=0)
        dQmJ[k] += _e.lmbda * intgrl[0]

    #########  HIGH
    for k from 1 <= k < pq_H:
        intgrl = _sciint.quad(intg1, _e.fH1, _e.fH2, args=(Fh, F2h, _e.IxH, _e.IIxH, _d.kH, k), epsrel=_e.epsrel, epsabs=0)
        dQmJ[pq_L+k] += _e.lmbda * intgrl[0]
    dQmJ *= -1
    dQmJ = dQmJ.reshape(pq_L + pq_H)

    return dQmJ

def mQlJ(Fnz, _d, _e):
    """
    Fnz    low and high AR + noizes
    """
    #  components 1 + p_L + 1 + p_H  (incl. noise term)
    pq_L = 1 + _d.kL;         pq_H = 1 + _d.kH
    Fnz= Fnz.reshape((pq_L + pq_H, 1))   #  need this to do outerproducts
    #  Fnz[0]    noize   Fnz[1]        ... Fnz[pq_L - 1]       dynamics
    #  Fnz[pq_L] noize   Fnz[pq_L + 1] ... Fnz[pq_L + pqH - 1] dynamics
    q2_L  = Fnz[0, 0];      q2_H  = Fnz[pq_L, 0]
    Fl  = Fnz[1:pq_L, 0];   Fh  = Fnz[pq_L+1:, 0]
    F2l = _N.dot(Fnz[1:pq_L, :], Fnz[1:pq_L, :].T)
    F2h = _N.dot(Fnz[pq_L + 1, :], Fnz[pq_L + 1, :].T)

    QlJ = -(_d.N / 2.)*(log(2*_N.pi*q2_L) + log(2*_N.pi*q2_H)) - (_e.CL - 2*_N.dot(Fl.T, _e.BL) + _N.dot(Fl.T, _N.dot(_e.AL, Fl))) / (2*q2_L) - (_e.CH - 2*_N.dot(Fh.T, _e.BH) + _N.dot(Fh.T, _N.dot(_e.AH, Fh))) / (2*q2_H)

    intgrl = _sciint.quad(intg0, _e.fL1, _e.fL2, args=(Fl, F2l, _e.IxL, _e.IIxL, _d.kL), epsrel=_e.epsrel, epsabs=0)
    QlJ  -= log(q2_L)*_e.lmbda * intgrl[0]

    #########  HIGH
    intgrl = _sciint.quad(intg0, _e.fH1, _e.fH2, args=(Fh, F2h, _e.IxH, _e.IIxH, _d.kH), epsrel=_e.epsrel, epsabs=0)
    QlJ  -= log(q2_H)*_e.lmbda * intgrl[0]
    QlJ *= -1
    return QlJ

def intg0(f, F0, FF0, Indx, IIndx, p):
    b = _b(f, F0, FF0, Indx, IIndx, p)
    if b == 0:
        print "!!!!!!!!!!!!!!!!!!!!!!   b == 0, intg0"

    return log(1. / b)

def intg1(f, F0, FF0, Indx, IIndx, p, k):
    b = _b(f, F0, FF0, Indx, IIndx, p)
    if b == 0:
        print "!!!!!!!!!!!!!!!!!!!!!!   b == 0, intg1"

    return _db(f, F0, Indx, p, k) / b

def intg2(f, F0, FF0, Indx, IIndx, p, k, l):
    b = _b(f, F0, FF0, Indx, IIndx, p)
    if b == 0:
        print "!!!!!!!!!!!!!!!!!!!!!!   b == 0, intg2"
    return (_db(f, F0, Indx, p, k) * _db(f, F0, Indx, p, l) -  2*cos(6.283185307179586*(k-l)*f) * b) / (b*b)
