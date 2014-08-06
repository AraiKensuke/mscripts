import scipy.stats as _ss
import scipy as _sci
import matplotlib.pyplot as _plt
import numpy as _N
import utilities as _U
import em as _em

def multiState3(cnts, blksz=50, start=None, stop=None, ICs=25):
    L = len(cnts)
    orig = _N.zeros(L, dtype=_N.int)
    if start == None:
        for i in range(L):
            if cnts[i] > 0:
                start = i
                break
    if stop == None:
        for i in range(L-1, -1, -1):
            if cnts[i] > 0:   
                stop = i   #  index of last trial with nonzero trial
                break
    ####   number of trials is (stop - start + 1)

    if (start == None) or (stop == None):
        return None, None, None, None, None
    roi  = _N.zeros(stop - start + 1, dtype=_N.int)
    orig = _N.zeros(stop - start + 1, dtype=_N.int)
    roi[:] = cnts[start:stop + 1]

    blks=(stop-start + 1)/blksz

    if blks <= 3:
        return None, None, None, None, None

    ssb     = _N.zeros((blks, 2), _N.int)  #  start stop blocks

    #  blks  rows x 2 cols
    abestfps =  _N.zeros((blks, 4))  #  first two center freqs, last two weights
    c2vs = _N.zeros(blks)
    ts   = _N.zeros(blks)
    lbsz = _N.log(blksz)
    for b in xrange(blks):
       i = b*blksz
       ssb[b, 0] = start + i
       ssb[b, 1] = start + i + blksz

       if _N.sum(roi[i:i + blksz]) == 0:   #  zero counts in all trials
          abestfps[b, 0:2] = _N.zeros(2)[:]
          abestfps[b, 2:4] = _N.array([1., 0.])[:]
          c2vs[b] = 1
       else:
          bestf, bestp, ll2, ll1 = _em.poissonEM_2comp(roi[i:i + blksz], ICs=ICs)
          
          BIC1 = -2*ll1 + lbsz
          BIC2 = -2*ll2 + 3*lbsz

          abestfps[b, 0:2] = bestf[:]
          abestfps[b, 2:4] = bestp[:]
          c2vs[b] = _N.std(roi[i:i + blksz])**2/_N.mean(roi[i:i + blksz])
          ts[b]   = 0
          if BIC2 < BIC1:
             ts[b]   = 1
    return abestfps, c2vs, ssb, ts, float(_N.sum(roi[:])) / (stop - start)
        #  Try two state fit test
