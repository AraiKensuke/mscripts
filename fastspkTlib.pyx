import numpy as _N

def spkTsInInterval(spkts, intvs, minIntv=500, maxIntv=None):
    """
    normally, spks in interval
    """
    bMaxIntvOK = False
    if maxIntv == None:
        bMaxIntvOK = True
    it0   = 0

    spktsA = []
    rtsA   = []
    L      = len(spkts)
    _spktsA= _N.zeros(10000, dtype=_N.int32)

    lenI   = len(intvs)
    for ic from 0 <= ic < lenI:
        intv = intvs[ic]
#    for intv in intvs:
        it = it0
        t0 = intv[0]
        t1 = intv[1]
        T  = t1 - t0

        if (T > minIntv) and (bMaxIntvOK or (T < maxIntv)):   #  interval itself is long
            spksInIntv = 0
            while it < L:
                if (spkts[it] >= t0) and (spkts[it] < t1):
                    _spktsA[spksInIntv] = spkts[it]
                    spksInIntv += 1
                if spkts[it] >= intv[1]:
                    break
                it += 1
            spktsA.append(_spktsA[0:spksInIntv].tolist())
            rtsA.append(spksInIntv / (float(T)) * 1000)
            it0 = it - 1
            if it0 < 0:
                it0 = 0
                
    return rtsA, spktsA
