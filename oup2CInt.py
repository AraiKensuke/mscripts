#  ctypes doc. 
#  http://docs.scipy.org/doc/numpy-1.3.x/user/c-info.python-as-glue.html

__all__ = ['OUP2']

import numpy as _N
import os
import ctypes as _ct
import numpy.ctypeslib as _ctl

#  _path = os.path.dirname('__file__')
#  _path = os.path.dirname(os.environ["LD_LIBRARY_PATH"])
#lib = _N.ctypeslib.load_library('oup2lib.dylib', _path)
if os.environ["OS_ARCH"] == "MacPPC":
    lib    = _ct.CDLL("oup2libPPC.dylib")    #  uses $LD_LIBRARY_PATH
else:
    lib    = _ct.CDLL("oup2libInt.dylib")    #  uses $LD_LIBRARY_PATH
#  read in C library file 


#  get/set attributes for C function called naisekiK
val = getattr(lib, 'OUP2')
val.restype = None   #  returns void
val.argtypes = [
    _ct.c_double, _ct.c_double, _ct.c_double,
    _ctl.ndpointer(_ct.c_double, ndim=1, 
                   flags='aligned,contiguous,writeable'), 
    _ctl.ndpointer(_ct.c_double, ndim=1, 
                   flags='aligned,contiguous,writeable'), 
    _ct.c_double]


#  wrapper function for naiseki.  Make sure variable types set correctly.  Note
#  that since function signatures are identical, we can use string libFuncName 
#  to specify actual function (naisekiK, naiseki2K or ...) to be called
def OUP2(tau, TEND, DT, x1, x2, p):
    requires = ['CONTIGUOUS', 'ALIGNED']
#  Return an ndarray of the provided type that satisfies requirements.
#  This function is useful to be sure that an array with the correct flags 
#  is returned for passing to compiled code (perhaps through ctypes)
    tau         = _N.require(tau, _ct.c_double)
    TEND        = _N.require(TEND, _ct.c_double)
    DT          = _N.require(DT, _ct.c_double)
    x1          = _N.require(x1, _ct.c_double, requires)
    x2          = _N.require(x2, _ct.c_double, requires)
    p           = _N.require(p, _ct.c_double)
    return lib.OUP2(tau, TEND, DT, x1, x2, p)
