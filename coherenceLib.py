import numpy as _N
import os
import ctypes as _ct
import numpy.ctypeslib as _ctl
from scipy.signal import hilbert



def smoothPhaseDiffEvo(ph1, ph2, maxPhaseDiffJump=1.0, diffLo=0.0, diffHi=6.28318530717959):
    requires = ['CONTIGUOUS', 'ALIGNED']
    L1 = _N.size(ph1)
    L2 = _N.size(ph2)

    if L1 > L2:
        ph1 = ph1[0:L2]
        L   = L2
        print "Input lengths not same.  Adjusting to length of 1st array"
    elif L2 > L1:
        ph2 = ph2[0:L1]
        L   = L1
        print "Input lengths not same.  Adjusting to length of 2nd array"
    else:
        L   = L1

    ph1   = _N.require(ph1, _ct.c_double, requires)
    ph2   = _N.require(ph2, _ct.c_double, requires)
    phaseDiff = _N.zeros(L, _ct.c_double)
    phaseDiff = _N.require(phaseDiff, _ct.c_double, requires)
    L   = _N.require(L, _ct.c_long)
    maxPhaseDiffJump   = _N.require(maxPhaseDiffJump, _ct.c_double)
    diffLo   = _N.require(diffLo, _ct.c_double)
    diffHi   = _N.require(diffHi, _ct.c_double)
    lib.smoothPhaseDiffEvo(ph1, ph2, phaseDiff, L, maxPhaseDiffJump, diffLo, diffHi)
    return phaseDiff


#  atan for 1 pair of (x, y)
def base_q4atan(x, y):
    if x == 0 and y > 0:
        atan =  0.5*_N.pi
    elif x == 0 and y < 0:
        atan = 1.5*_N.pi
    elif y >= 0 and x > 0:#Q1
        atan = _N.arctan(y / x)
    elif y > 0 and x <= 0:#Q2
        atan = _N.arctan(y / x) + _N.pi
    elif y <= 0 and x < 0:#Q3
        atan = _N.arctan(y / x) + _N.pi
    elif y < 0 and x >= 0:
        atan = _N.arctan(y / x) + 2*_N.pi
    else:
        atan = -100
    return atan

#  define function "atan" to be an element-by-element operation
atan = _N.vectorize(base_q4atan)

#  to a time series of phase
def toAmpPhase(x):
    xH = hilbert(x)
    return _N.real(xH)**2 + _N.imag(xH)**2, atan(_N.real(xH), _N.imag(xH))

def avgPhaseAndOP(phaseDiffs):
    N = _N.size(phaseDiffs)
    rePt = 0
    imPt = 0

    for n in range(0, N):
        rePt += _N.cos(phaseDiffs[n])
        imPt += _N.sin(phaseDiffs[n])
        
    RR = (1.0/N) * _N.sqrt(rePt*rePt + imPt*imPt)
    TTHETA = _N.arctan(imPt / rePt);  #  in QUAD IV, I

    if rePt < 0 and imPt > 0:
        TTHETA += _N.pi
    elif imPt < 0 and rePt < 0: 
        TTHETA += _N.pi
    elif rePt > 0 and imPt < 0:
        TTHETA += 2*_N.pi

    return RR, TTHETA
    
