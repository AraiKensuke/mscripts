import numpy as _N
import os
import ctypes as _ct
import numpy.ctypeslib as _ctl
from scipy.signal import hilbert

# #  Input N row by 2 col of data.  Correlation == 1 is simple to interpret.
# #  However, correlation ~ 0 is more difficult.  Is it because col. 1 and
# #  col. 2 are statistically independent, or is it because there is some
# #  phase delay between the two signals?  slideCorr adds a range of delays
# #  btwn. col. 1 and col. 2, and returns the delay that resulted in the 
# #  largest correlation value as well as this value of correlation
# #  
# #  col.1 :  RRRRRRRRRRRRRRRRRR
# #
# #  col. 2:  mmmmOOOOOOOOOOOOOO
# #  col. 2': OOOOOOOOOOOOOOmmmm    (delay < 0 case)
# #  
# #  for delay < 0, we see that earlier col. 1 events are now being
# #  compared to later col. 2 events.  ==>  col. 1 delayed wrt. col. 2
# def slideCorr(twoColDat, pmNSteps):
#     temp = _N.zeros(_N.shape(twoColDat))
#     temp[:,0] = twoColDat[:,0]
#     dataLen = _N.shape(temp)[0]
#     maxCorr = -2
#     maxCorrDel = -pmNSteps

#     for delay in range(-pmNSteps, pmNSteps + 1):
#         if delay < 0:
#             mdelay = -delay
#             temp[mdelay:dataLen,1] = twoColDat[0:dataLen - mdelay, 1]
#             temp[0:mdelay,1] = twoColDat[dataLen - mdelay:dataLen, 1]
#         elif delay >= 0:
#             temp[0:dataLen - delay,1] = twoColDat[delay:dataLen, 1]
#             temp[dataLen - delay:dataLen,1] = twoColDat[0:delay, 1]

#         corr = _N.corrcoef(temp.T)[0, 1]
#         if corr > maxCorr:
#             maxCorr = corr
#             maxCorrDel = delay

#     return maxCorr, maxCorrDel


# #  ctypes doc. 
# #  http://docs.scipy.org/doc/numpy-1.3.x/user/c-info.python-as-glue.html

# __all__ = ['correntropy', 'smoothPhaseDiffEvo']

# #_path = os.path.dirname('__file__')
# #_path = os.path.dirname(os.environ["LD_LIBRARY_PATH"])
# #lib = _N.ctypeslib.load_library('correntropyCcode', _path)
# #  read in C library file 
# if os.environ["OS_ARCH"] == "MacPPC":
#     lib    = _ct.CDLL("coherenceCcodePPC.dylib")    #  uses $LD_LIBRARY_PATH
# else:
#     lib    = _ct.CDLL("coherenceCcodeInt.dylib")    #  uses $LD_LIBRARY_PATH


# #  get/set attributes for C function called naisekiK
# val = getattr(lib, 'correntropy')
# val.restype = _ct.c_double   #  returns void
# val.argtypes = [_ctl.ndpointer(_ct.c_double, ndim=1, 
#                                flags='aligned,contiguous'), 
#                 _ctl.ndpointer(_ct.c_double, ndim=1, 
#                                flags='aligned,contiguous'), 
#                 _ct.c_int, _ct.c_double]


# #  wrapper function for naiseki.  Make sure variable types set correctly.  Note
# #  that since function signatures are identical, we can use string libFuncName 
# #  to specify actual function (naisekiK, naiseki2K or ...) to be called
# def correntropy(x, y, sig2):
#     requires = ['CONTIGUOUS', 'ALIGNED']
# #  Return an ndarray of the provided type that satisfies requirements.
# #  This function is useful to be sure that an array with the correct flags 
# #  is returned for passing to compiled code (perhaps through ctypes)
#     x           = _N.require(x, _ct.c_double, requires)
#     y           = _N.require(y, _ct.c_double, requires)
#     sig2        = _N.require(sig2, _ct.c_double)
#     return lib.correntropy(x, y, _N.size(x), sig2)

# #  get/set attributes for C function called naisekiK
# val = getattr(lib, 'smoothPhaseDiffEvo')
# val.restype = None   #  returns void
# val.argtypes = [_ctl.ndpointer(_ct.c_double, ndim=1, 
#                                flags='aligned,contiguous'), 
#                 _ctl.ndpointer(_ct.c_double, ndim=1, 
#                                flags='aligned,contiguous'), 
#                 _ctl.ndpointer(_ct.c_double, ndim=1, 
#                                flags='aligned,contiguous'), 
#                 _ct.c_long, _ct.c_double, _ct.c_double, _ct.c_double]

#  from two time evolutions of phase, take the difference.  however, 
#  if ph1 = (6.0, 6.2, 0.2, 0.4), ph2 = (5.9, 6.0, 6.2, 0.2)
#  diff = (0.1, 0.2, -6, 0.2)  NOT GOOD, want  (0.1, 0.2, -6+2pi=0.28, 0.2)
#  use assumption that evolution of phase diff is small
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
  if y >= 0 and x > 0:#Q1
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
def toPhase(x):
    xH = hilbert(x)
    return atan(_N.real(xH), _N.imag(xH))


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
    
