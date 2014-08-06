import scipy.stats as _ss
import scipy as _sci
import matplotlib.pyplot as _plt
import numpy as _N

def poi_pdf(l, k):
    return (l**k * _N.exp(-l)) / _sci.factorial(k)

def poi_ch2test(cts):
    _plt.ioff()
    TRIALS    = len(cts)

    rareLimL = 0
    rareLimH = 0
    obsLam   = _N.mean(cts)

    i = int(obsLam)
    while True:
        if poi_pdf(obsLam, i)*TRIALS < 1:
            rareLimH = i - 1
            break   #  inclusive of rareLim and up
        i += 1

    i = int(obsLam)
    while True and (i >= 0):
        if poi_pdf(obsLam, i)*TRIALS < 1:
            rareLimL = i + 1
            break   #  inclusive of rareLim and up
        i -= 1

    #  no bin has < 1 expcted events

    expctd   = _N.zeros(rareLimH + 1)
    for n in range(rareLimH):
        expctd[n] = poi_pdf(obsLam, n)*TRIALS
    expctd[rareLimH] = TRIALS - _N.sum(expctd[0:rareLimH])

    #  poipdf[rareLim] is at < 1.  This will be last square.
    #  so this is rareLim + 1 objects
    maxInd              = max(cts)  #  0 based index
    vals,  bins,  objs  = _plt.hist(cts, bins=range(0, maxInd + 2, 1))
    #  # of categories:  rareLim + 1    0..4 -> 1..5  (5 cats).  
    #  k == # of classes

    #  shortened version

    #  if rareLim is the last spot, then length is rareLim + 1
    #  highest count in cts is len(vals) - 1
    #  maxLim == len(vals) - 1   
    #  to accomodate last index rareLim, we need size rareLim + 1
    #  
    svals                =  _N.zeros(rareLimH + 1)
    if maxInd <= rareLimH:
        svals[0:maxInd + 1] =  vals[0:maxInd + 1]
        
    else:
        svals[0:rareLimH] =  vals[0:rareLimH]
        svals[rareLimH]       =  _N.sum(vals[rareLimH:])

    svals[rareLimL]      = _N.sum(vals[0:rareLimL+1])
    expctd[rareLimL]     = _N.sum(expctd[0:rareLimL+1])

    k                   = rareLimH - rareLimL + 1  # index of last element if index from 1
    # [0, 1], [1, 2], ... [k-2, k-1]
    #   # classes (counts) [0, 1, 2, k-2]   (k - 1) classes
    #  1...k inclusive is k classes

    chi2   = 0

    for i in xrange(rareLimL, rareLimH + 1):
        o    =  svals[i]   
        e    =  expctd[i]
        chi2  +=  (o-e)**2/e
        i    += 1

    edf  = k - 2
    pv = 1 - _ss.chi2.cdf(chi2, edf)

    _plt.ion()

    return pv
