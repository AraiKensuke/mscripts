import numpy as _N

def tf():
    cdef double j = 0
    cdef int i
    for i from 0 <= i < 1000000:
        j = _N.sin(i)

    print j


