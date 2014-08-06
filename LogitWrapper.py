import _LogitWrapper as __lw
import numpy as _N

## Draw PG(h, z)
##------------------------------------------------------------------------------
#  rpg.gamma <- function(num=1, h=1, z=0.0, trunc=200)

## Draw PG(n, z) where n is a natural number.
##------------------------------------------------------------------------------
#  rpg.devroye <- function(num=1, n=1, z=0.0)

## Draw PG(h, z) where h is >= 1.
##------------------------------------------------------------------------------
#  rpg.alt <- function(num=1, h=1, z=0.0)

## Draw PG(h, z) using SP approx where h is >= 1.
##------------------------------------------------------------------------------
#  rpg.sp <- function(num=1, h=1, z=0.0, track.iter=FALSE)

def rpg_devroye(n, z, num=1, out=None):
    arr_n = (type(n) == _N.ndarray) or (type(n) == list) or (type(n) == tuple)
    arr_z = (type(z) == _N.ndarray) or (type(z) == list) or (type(z) == tuple)

    if (num > 1) or arr_n or arr_z:
        if out == None:
            out = _N.empty(num, dtype=_N.float)
        if not arr_n:
            int_n = (type(n) == _N.int16) or (type(n) == _N.int32) or (type(n) == _N.int64) or (type(n) == int)
            if int_n:
                n = _N.ones(num, dtype=_N.int32) * n
            else:
                #  throw an exception
                pass
        if not arr_z:
            z = _N.ones(num, dtype=_N.float) * z
        __lw.rpg_devroyeM(out, n, z, num)
        return out
    else:  #  both n, z not arrays
        int_n = (type(n) == _N.int16) or (type(n) == _N.int32) or (type(n) == _N.int64) or (type(n) == int)
        if not int_n:
            # throw an exception
            pass
        elif n < 0:
            # throw an exception
            pass
        if z < 0:
            # throw an exception
            pass
        z = 1.0 * z  #  force to float
        x = __lw.rpg_devroyeS(n, z)
        return x

def rpg_gamma(h, z, num=1, trunc=200):
    #  no num.  
    arr_h = (type(h) == _N.ndarray) or (type(h) == list) or (type(h) == tuple)
    arr_z = (type(z) == _N.ndarray) or (type(z) == list) or (type(z) == tuple)
    lh    = 1
    lz    = 1
    if arr_h:
        lh = len(h)
    if arr_z:
        lz = len(z)

    if (num > 1) or arr_h or arr_z:
        if (num == 1) and (lh > 1):
            num = lh
        if (num == 1) and (lz > 1):
            num = lz
        if (lh > 1) and (lz > 1) and (lh != lz):
            #  throw an exception
            pass

        retArr = _N.empty(num, dtype=_N.float)
        if not arr_h:
            float_h = (type(h) == _N.float16) or (type(h) == _N.float32) or (type(h) == _N.float64) or (type(h) == float)
            if float_h:
                h = _N.ones(num, dtype=_N.float) * h
            else:
                h = _N.ones(num, dtype=_N.float) * h
        if not arr_z:
            z = _N.ones(num, dtype=_N.float) * z
        __lw.rpg_gamma(retArr, h, z, num, 200)
        return retArr
    else:  #  both n, z not arrays
        int_h = (type(h) == _N.int16) or (type(h) == _N.int32) or (type(h) == _N.int64) or (type(h) == int)
        if not int_h:
            # throw an exception
            pass
        elif h < 0:
            # throw an exception
            pass
        if z < 0:
            # throw an exception
            pass
        z = 1.0 * z  #  force to float
        x = __lw.rpg_gamma(h, z)
        return x

def rpg_alt(h, z, num=1, trunc=200):
    arr_h = (type(h) == _N.ndarray) or (type(h) == list) or (type(h) == tuple)
    arr_z = (type(z) == _N.ndarray) or (type(z) == list) or (type(z) == tuple)

    if (num > 1) or arr_h or arr_z:
        retArr = _N.empty(num, dtype=_N.float)
        if not arr_h:
            int_h = (type(h) == _N.int16) or (type(h) == _N.int32) or (type(h) == _N.int64) or (type(h) == int)
            if int_h:
                h = _N.ones(num, dtype=_N.int32) * h
            else:
                #  throw an exception
                pass
        if not arr_z:
            z = _N.ones(num, dtype=_N.float) * z
        __lw.rpg_alt(retArr, h, z, num)
        return retArr
    else:  #  both n, z not arrays
        int_h = (type(h) == _N.int16) or (type(h) == _N.int32) or (type(h) == _N.int64) or (type(h) == int)
        if not int_h:
            # throw an exception
            pass
        elif h < 0:
            # throw an exception
            pass
        if z < 0:
            # throw an exception
            pass
        z = 1.0 * z  #  force to float
        x = __lw.rpg_alt(h, z)
        return x
