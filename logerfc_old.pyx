import pickle as _pkl
import numpy as _N
import scipy.stats as _ss

cdef extern from "math.h":
    double exp(double)
    double sqrt(double)
    double log(double)
    double abs(double)

_logerfc = None
_x0      = 0
_x1      = None
_R       = None
_dx      = None
#  useful constants
_sqrt2   = sqrt(2)
_twpi    = 2*_N.pi

def init():
    global _logerfc, _x1, _R, _dx

    if _logerfc == None:
        fpk = open("/Users/arai/usb/nctc/mscripts/logerfc.dat", "r")
        dat = _pkl.load(fpk)
        _logerfc = dat
        _R, c     = _logerfc.shape
        _x1      = _logerfc[_R-1, 0]
        _dx      = _logerfc[1, 0] - _logerfc[0, 0]
        fpk.close()

def at(x):
    global _x1, _logerfc, _R, _dx
    abx = abs(x)
    
    if abx > _x1:
        #print "abx  %(abx).3e   >  _x1  %(x1).3e" % {"abx" : abx, "x1" : _x1}
        return None
        
    n   = int(abx / _dx)
    rx  = (abx - _logerfc[n, 0]) / _dx

    if n < _R - 1:
        lval = _logerfc[n, 1] + rx * (_logerfc[n+1, 1] - _logerfc[n, 1])
    else:
        lval = _logerfc[n, 1]

    if x < 0:
        #  log(2 - val) = log(val x [2/val - 1])
        lval = log(2 - exp(lval))

    return lval

def trncNrmNrmlz(a, b, u, sg):
    #   single variable Normal distribution
    #  Value of density is 
    #  _N.exp(-(xs - u)**2 / (2*sg2) - lNC)

    #  Normalization is 
    global _sqrt2, _twpi
    anmz = (a-u) / (_sqrt2*sg)
    bnmz = (b-u) / (_sqrt2*sg)
    sg2  = sg*sg

    #################################
    if (anmz < 0) and (bnmz <= 0):
        tmp  = anmz
        anmz = -bnmz
        bnmz = -tmp
        #print "A   anmz  %(1).3e    bnmz  %(2).3e" % {"1" : anmz, "2" : bnmz}
        eps = exp(at(bnmz) - at(anmz))
        lNC = log(0.5) + log(sqrt(_twpi*sg2)) + at(anmz) - eps - 0.5*eps*eps - 0.166666666*eps*eps*eps - 0.04166666666*eps*eps*eps*eps
    elif (anmz >= 0) and (bnmz > 0):
        #print "B  anmz  %(1).3e    bnmz  %(2).3e" % {"1" : anmz, "2" : bnmz}
        eps = exp(at(bnmz) - at(anmz))
        lNC = log(0.5) + log(sqrt(_twpi*sg2)) + at(anmz) - eps - 0.5*eps*eps - 0.166666666*eps*eps*eps - 0.04166666666*eps*eps*eps*eps
    #################################
    elif (anmz < 0) and (bnmz > 0):
        NC = sqrt(2*_N.pi*sg2) * (_ss.norm.cdf((b - u) / sg) - \
                                     _ss.norm.cdf((a - u) / sg))
        lNC= log(NC)
    else:
        print "woops"
        print anmz
        print bnmz

    return lNC
