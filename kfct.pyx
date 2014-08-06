import numpy as _N
cimport numpy as _N   #  need to add include path  /Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/numpy/core/include
#import tokyo as _tky
#cimport tokyo as _tky

"""
c functions
"""
cdef extern from "math.h":
    double exp(double)
    double sqrt(double)
    double log(double)
    double abs(double)

#  estimate mode of posterior
cimport cython
@cython.boundscheck(False)
def f_and_df(x, *args):
    #    args=(F, flt_x[:, 0, n-1], _N.linalg.inv(pr_V[:, :, n]), beta, u, dN, n)
    #  can't pass x as column vector.  "ValueError: ... too deep"
#    print x
    k   = args[0]
    xnnm= args[1]
    Vnnm= args[2]
    Vinnm= args[3]
    u    = args[4]
    y    = args[5]
    r    = args[6]
    x    = x.reshape(k, 1)

    ef    = exp(u + x[0, 0])   # exp using filtered state
    pspk  = ef / (1 + ef)
    #  xnmnm is formatted so xnmnm[0] is most recent value
#    f  = -_tky.dgemm(Vinnm, x - xnnm) + beta * (dN - pspk)
    f  = -_N.dot(Vinnm, x - xnnm) + (y - (r + y) * pspk)
    df = -Vinnm - (r + y) * (ef / ((1+ef)*(1+ef)))
    f = f.reshape(k)

    return f, df
