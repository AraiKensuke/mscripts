import pickle as _pkl
import numpy as _N
import scipy.stats as _ss

_ar2propArea = None
_x0      = 0
_x1      = None
_R       = None
_dx      = None
#  useful constants
_sqrt2   = _N.sqrt(2)
_twpi    = 2*_N.pi

def init(sg):
    global _ar2propArea, _x1, _R, _dx

    if _ar2propArea == None:
        fpk = open(prcmpFN("ar2IntgArea,%.3f" % sg), "r")
        dat = _pkl.load(fpk)
        _ar2propArea = dat
        _R, c     = _ar2propArea.shape
        _x1      = _ar2propArea[_R-1, 0]
        _dx      = _ar2propArea[1, 0] - _ar2propArea[0, 0]
        fpk.close()

def at(x):
    global _x1, _ar2propArea, _R, _dx
    abx = _N.abs(x)
    
    if abx > _x1:
        return None
        
    n   = int(abx / _dx)
    rx  = (abx - _ar2propArea[n, 0]) / _dx

    if n < _R - 1:
        lval = _ar2propArea[n, 1] + rx * (_ar2propArea[n+1, 1] - _ar2propArea[n, 1])
    else:
        lval = _ar2propArea[n, 1]

    if x < 0:
        #  log(2 - val) = log(val x [2/val - 1])
        lval = _N.log(2 - _N.exp(lval))

    return lval
