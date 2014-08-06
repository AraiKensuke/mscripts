import scipy.stats as _ss
import scipy as _sci
import matplotlib.pyplot as _plt
import numpy as _N
import utilities as _U

def poi_pdf(l, karr):
   rv = _ss.poisson(l)
   return rv.pmf(karr)

def poi_ch2test(cts):
    if _N.sum(cts) == 0:  #  all of them 0
        return 0.5    #  Poisson with rate 0
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
    nbins = maxInd - 0 + 1

    hf, low, bs, out = _ss.histogram(cts, numbins=nbins, defaultlimits=(0, maxInd + 1))

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
        svals[0:maxInd + 1] =  hf[0:maxInd + 1]
        
    else:
        svals[0:rareLimH] =  hf[0:rareLimH]
        svals[rareLimH]       =  _N.sum(hf[rareLimH:])

    svals[rareLimL]      = _N.sum(hf[0:rareLimL+1])
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

    return pv

#  cumulative fraction
def cumfrac(x, staircase=False, countdat=False, histogram=False, bins=None, binsAlignLeft=False):
    """
    cumulative fraction
    types of data:

    continuous data        sort data, then assign size rank to sorted points
    countdata              make a histogram of data first
    histogram              cumulatively add values of bins

    cum frac looks like a staircase when plotted w/ lines
    staircase=True will include the 
    cnts   if data is a list of counts (like spks per trial), we will most likely have many instances of, ie the number 4, in our data.  We just need to consider all 

    sanity check
    cnts = [1, 2, 3, 4, 5]  or [1, 2, 3, 4, 5, 3]
    cf = _ks.cumfrac(cnts, countdat=True, staircase=True)    
    plot(cf[:, 0], cf[:, 1])
    """

    if countdat:
        lo = min(x)
        hi = max(x)
        nbins = hi - lo + 1
        hf, low, bs, out = _ss.histogram(x, numbins=nbins, defaultlimits=(lo, hi + 1))

        datn = int(_N.sum(hf))
        if not staircase:
            cf   = _N.zeros((hi - lo + 1, 2))
            tot  = 0
            for i in xrange(len(hf)):
                tot += hf[i]
                cf[i, 0] = lo + i
                cf[i, 1] = float(tot) / datn
        else:
            cf   = _N.zeros(((hi - lo + 1)*2, 2))
            tot  = 0
            for i in xrange(len(hf)):
                cf[2*i, 0]     = lo + i
                cf[2*i + 1, 0] = lo + i
                cf[2*i, 1]     = float(tot) / datn
                tot            += hf[i]
                cf[2*i + 1, 1] = float(tot) / datn

    else:
        sx = _N.sort(x)
        N  = len(sx)

        if not staircase:  
            cf = _N.zeros((N, 2))
            cf[:, 0] = sx[:]
            cf[:, 1] = _N.linspace(0, 1, N, endpoint=False) + 1./N
        else:
            cf = _N.zeros((2*N, 2))

            yvals = _N.linspace(0, 1, N, endpoint=False)
            for n in xrange(N):
                cf[2*n, 0]      =  sx[n]
                cf[2*n + 1, 0]  =  sx[n]
                cf[2*n, 1]      =  yvals[n]
                cf[2*n + 1, 1]  =  yvals[n] + 1./N
                
    if histogram:
        tot = _N.sum(x)
        N   = len(x)
        cf  = _N.zeros((N, 2))
        ccf = 0
        if bins == None:
            bins = range(N + 1)
        for b in xrange(N):
            cf[b, 0] = bins[b + 1]
            cf[b, 1] = ccf
            ccf += x[b]
        cf[N - 1, 0] = bins[b + 1]
        cf[N - 1, 1] = tot

    return cf

def percentile(x, countdat=False):
    if countdat:
        lo = min(x)
        hi = max(x)
        nbins = hi - lo + 1
        hf, low, bs, out = _ss.histogram(x, numbins=nbins, defaultlimits=(lo, hi + 1))

        datn = int(_N.sum(hf))
        pctl   = _N.zeros((hi - lo + 1, 2))
        tot  = 0
        for i in xrange(len(hf)):
            tot += hf[i]
            pctl[i, 0] = lo + i
            pctl[i, 1] = float(tot + 1) / float(datn + 1)

    else:
        sx = _N.sort(x)
        N  = len(sx)

        pctl = _N.zeros((N, 2))

        for n in xrange(N):
            pctl[n, 0]      =  sx[n]
            pctl[n, 1]      =  float(n + 1) / float(N + 1)
            
    return pctl

def subtractDiffSupport(xy1, xy2):
    """
    subtract two functions that have a different support.
    occurs for example when the functions are the result of random sampling
    """

    p2    = 0
    ix    = 0
    lastxy1_y = 0
    lastxy2_y = 0
    
    L1 = len(xy1[:, 0])
    L2 = len(xy2[:, 0])
    subt_x = []
    subt_y = []
    for ix in xrange(L1):
        xy1_x = xy1[ix, 0]
        xy1_y = xy1[ix, 1]

        while (p2 < L2) and (xy2[p2, 0] < xy1_x):
            subt_x.append(xy2[p2, 0])
            subt_y.append(lastxy1_y - xy2[p2, 1])
            lastxy2_y     = xy2[p2, 1]
            p2 += 1

        subt_x.append(xy1_x)
        subt_y.append(xy1_y - lastxy2_y)
        lastxy1_y = xy1_y

    subtracted = _N.array([subt_x, subt_y]).T

    return subtracted

def subtract_discrete_cf(xy1, xy2):
    """
    subtract two functions that have (almost) same support.  
    """

    lo1 = int(min(xy1[:, 0]))
    hi1 = int(max(xy1[:, 0]))
    lo2 = int(min(xy2[:, 0]))
    hi2 = int(max(xy2[:, 0]))

    LO  = min([lo1, lo2])
    HI  = max([hi1, hi2])

    y1  = _N.zeros(HI - LO + 1)
    y2  = _N.zeros(HI - LO + 1)
    x1  = _N.zeros(HI - LO + 1)
    x2  = _N.zeros(HI - LO + 1)


    y1[0:lo1-LO] = 0
    
    y1[lo1-LO:lo1-LO+(hi1-lo1+1)] = xy1[:, 1]
    y1[hi1+1-LO:HI-LO+1] = 1

    y2[0:lo2-LO] = 0
    y2[lo2-LO:lo2-LO+(hi2-lo2+1)] = xy2[:, 1]
    y2[hi2+1-LO:HI-LO+1] = 1

    return range(LO, HI + 1), y1 - y2


def poissonEM_2comp(cts, ICs=20):
    """
    fit data cts with mixture Poisson distribution using EM
    EM depends on initial conditions.  ICs is # of initial conditions to try
    """
    __inf__ = float("inf")
    N       = len(cts)
    K       = 2
    pik     = _N.zeros(K)
    lamk    = _N.zeros(K)
    bestf   = _N.zeros(K)
    bestp   = _N.zeros(K)

    K       = 2
    pik     = _N.zeros(K)
    lamk    = _N.zeros(K)
    bestf   = _N.zeros(K)
    bestp   = _N.zeros(K)
    ICs     = 25

    maxLL = -10000000.

    for ic in xrange(ICs):
        pik[0]  = _N.random.rand()
        pik[1]  = 1 - pik[1]

        lamk[0] = _N.mean(cts) * (0.7 + 0.6*_N.random.rand())
        lamk[1] = _N.mean(cts) * (0.7 + 0.6*_N.random.rand())

        Nk      = _N.zeros(K)
        resp    = _N.zeros((N, K))

        llo     = -10000000000.

        while True:
            #  build responsibility array
            for n in xrange(N):
                bot = 0

                for k in xrange(K):
                    bot += pik[k] * poi_pdf(lamk[k], cts[n])

                for k in xrange(K):
                    resp[n, k] = (pik[k] * poi_pdf(lamk[k], cts[n])) / bot
            #  build responsibility array        
            for k in xrange(K):
                Nk[k] = 0
                for n in xrange(N):
                    Nk[k] += resp[n, k]

            #  update params
            for k in xrange(K):
                lamk[k] = 0
                for n in xrange(N):
                    lamk[k] += cts[n] * resp[n, k]
                lamk[k] /= Nk[k]

                pik[k] = Nk[k] / float(N)

            #  log likelihood
            ll   = 0
            for n in xrange(N):
                logarg = 0
                for k in xrange(K):
                    logarg += pik[k] * poi_pdf(lamk[k], cts[n])

                ll += _N.log(logarg)

            if (_N.abs(llo - ll) / _N.abs(_N.mean(llo + ll))) < 0.001:
                break
            llo = ll

            if (ll != ll) or (ll == __inf__):
                #  reinitialize, hit floating point error
                print "hit float error, reinit params, kstat.poissonEM_2comp()"
                pik[0]  = _N.random.rand()
                pik[1]  = 1 - pik[1]
                lamk[0] = _N.mean(cts) * (0.7 + 0.6*_N.random.rand())
                lamk[1] = _N.mean(cts) * (0.7 + 0.6*_N.random.rand())

            if ll > maxLL:
                maxLL = ll
                bestf[0] = lamk[0]
                bestf[1] = lamk[1]
                bestp[0] = pik[0]
                bestp[1] = pik[1]
        
        #  I want lamk[0] to be the lower one
        if bestf[0] > bestf[1]:
            temp    = bestf[0]
            bestf[0] = bestf[1]
            bestf[1] = temp
            temp    = bestp[0]
            bestp[0]  = bestp[1]
            bestp[1]  = temp

    return bestf, bestp

def rmOutliers(arr):
    return arr
#    return srtd[1:len(arr)-1]

def multiState(cnts, blksz=50, winsz=5, start=None, stop=None):
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
                stop = i
                break

    if (start == None) or (stop == None):
        return None, None, None, None, None, None
    roi  = _N.zeros(stop - start, dtype=_N.int)
    orig = _N.zeros(stop - start, dtype=_N.int)
    roi[:] = cnts[start:stop]

    hlfblksz = blksz/2
    blks=(stop-start)/hlfblksz - 1
    smpv_f = 0
    smpv_b = 0
    cvs_f  = []
    cvs_b  = []
    tr_f   = 0
    tr_b   = 0

    if blks <= 4:
        return None, None, None, None, None, None

    ssb     = _N.zeros((blks, 4), _N.int)  #  start stop blocks
    trdb    = [[], []]
    nonzBlksFW= blks
    nonzBlksBW= blks
    for b in xrange(blks):
        fwi = b*hlfblksz
        bwi = (stop-start)-(b+2)*hlfblksz
        rmof= roi[fwi:fwi+blksz]
        rmob= roi[bwi:bwi+blksz]

        chi2pvwoo_f = poi_ch2test(rmof)
        chi2pvwoo_b = poi_ch2test(rmob)

        ssb[b, 0] = start + fwi
        ssb[b, 1] = start + fwi + blksz
        ssb[b, 2] = start + bwi
        ssb[b, 3] = start + bwi + blksz

        if chi2pvwoo_f < 0.01:
            smpv_f += 1
        if _N.sum(rmof) == 0:   #  all 0
            nonzBlksFW -= 1
        else:
            if trend(roi[fwi:fwi+blksz], winsz=winsz):
                trdb[0].append(b)
                if chi2pvwoo_f:
                    tr_f += 1

        if chi2pvwoo_b < 0.01:
            smpv_b += 1
        if _N.sum(rmob) == 0:   #  all 0
            nonzBlksBW -= 1
        else:
            if trend(roi[bwi:bwi+blksz], winsz=winsz):
                trdb[1].append(b)
                if chi2pvwoo_b:
                    tr_b += 1

    #  blks  rows x 4 cols
    abestfps =  _N.zeros((blks, 4))
#    if (smpv_f + smpv_b > nonzBlksFW - smpv_f + nonzBlksBW - smpv_b):
    c2vs = _N.zeros(blks)
    for b in xrange(blks):
       i = b*hlfblksz
       if _N.sum(roi[i:i + blksz]) == 0:
          abestfps[b, 0:2] = _N.zeros(2)[:]
          abestfps[b, 2:4] = _N.array([1., 0.])[:]
          c2vs[b] = 1
       else:
          bestf, bestp = poissonEM_2comp(roi[i:i + blksz])
          abestfps[b, 0:2] = bestf[:]
          abestfps[b, 2:4] = bestp[:]
          c2vs[b] = _N.std(roi[i:i + blksz])**2/_N.mean(roi[i:i + blksz])

    return abestfps, c2vs, smpv_f + smpv_b, nonzBlksFW - smpv_f + nonzBlksBW - smpv_b, ssb, trdb
        #  Try two state fit test

def multiState2(cnts, blksz=50, start=None, stop=None):
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
                stop = i
                break

    if (start == None) or (stop == None):
        return None, None, None, None
    roi  = _N.zeros(stop - start, dtype=_N.int)
    orig = _N.zeros(stop - start, dtype=_N.int)
    roi[:] = cnts[start:stop]

    hlfblksz = blksz/2
    blks=(stop-start)/hlfblksz - 1

    if blks <= 4:
        return None, None, None, None

    ssb     = _N.zeros((blks, 4), _N.int)  #  start stop blocks

    #  blks  rows x 4 cols
    abestfps =  _N.zeros((blks, 4))
    c2vs = _N.zeros(blks)
    for b in xrange(blks):
       fwi = b*hlfblksz
       bwi = (stop-start)-(b+2)*hlfblksz
       i = fwi
       ssb[b, 0] = start + fwi
       ssb[b, 1] = start + fwi + blksz
       ssb[b, 2] = start + bwi
       ssb[b, 3] = start + bwi + blksz

       if _N.sum(roi[i:i + blksz]) == 0:
          abestfps[b, 0:2] = _N.zeros(2)[:]
          abestfps[b, 2:4] = _N.array([1., 0.])[:]
          c2vs[b] = 1
       else:
          bestf, bestp = poissonEM_2comp(roi[i:i + blksz])
          abestfps[b, 0:2] = bestf[:]
          abestfps[b, 2:4] = bestp[:]
          c2vs[b] = _N.std(roi[i:i + blksz])**2/_N.mean(roi[i:i + blksz])
    ts = twostate(abestfps, 0.05)

    return abestfps, c2vs, ssb, ts
        #  Try two state fit test

def c2v(cnts, blksz=50, start=None, stop=None):
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
                stop = i
                break

    if (start == None) or (stop == None):
        return None, None, None, None, None, None
    roi  = _N.zeros(stop - start, dtype=_N.int)
    orig = _N.zeros(stop - start, dtype=_N.int)
    roi[:] = cnts[start:stop]

    hlfblksz = blksz/2
    blks=(stop-start)/hlfblksz - 1

    #  blks  rows x 4 cols
    c2vs = _N.zeros(blks)
    for b in xrange(blks):
       i = b*hlfblksz
       if _N.sum(roi[i:i + blksz]) == 0:       #  not firing in any of these trials
          c2vs[b] = 0
       else:
          c2vs[b] = _N.std(roi[i:i + blksz])**2/_N.mean(roi[i:i + blksz])

    return _N.mean(c2vs)

def trend(ct, winsz=5, useBott=2):
    """
    does data ct have a trend?  look at 
    useBott == use lowest "useBott" numbers per block to assess trend
    """
    blks  = len(ct)/winsz
    xs    = range(blks)
    emins = []

    if useBott > 1:
        for i in xrange(blks):
            srtd = _N.sort(ct[winsz*i:winsz*(i+1)])
            emins.append(_N.sum(srtd[0:useBott]))
    else:
        for i in xrange(blks):
            emins.append(min(ct[winsz*i:winsz*(i+1)]))   
    (a_s,b_s,r,tt,stderr)=_ss.linregress(xs, emins)

    eas = abs(a_s)

    #  if slope is 0, don't even do test, just return False.
    #  often helps when count rate is very low, and most are just 0
    if eas == .0:
        return False

    TRIALS = 100
    bgr    = 0   #  slope of data surrogate bigger than data
    smins  = _N.zeros(blks)
    for tr in xrange(TRIALS):
        #  scramble
        sct = _U.shuffle(ct)

        #  sort is by default ascending
        if useBott > 1:
            for i in xrange(blks):
                srtd = _N.sort(sct[winsz*i:winsz*(i+1)])
                smins[i]  = _N.sum(srtd[0:useBott])
        else:
            for i in xrange(blks):
                smins[i] = (min(sct[winsz*i:winsz*(i+1)]))
        (a_s,b_s,r,tt,stderr)=_ss.linregress(xs, smins)

        ma_s = abs(a_s)

        if (eas < ma_s):
            bgr += 1
        
    if float(bgr) / TRIALS < 0.05:
        return True
    return False

def threshold(cnts, abestfps, blksz=50):
   """
   turn 2-state spk count data into a 0, 1 sequence.  
   cnts len L, return sequence hilos and thresh integer multiples of blksz
   abestfps is overlapping blocks of size blksz, slid every hlfblksz 
   """
   L = len(cnts)
   orig = _N.zeros(L, dtype=_N.int)
   for i in range(L):
      if cnts[i] > 0:
         start = i
         break
   for i in range(L-1, -1, -1):
      if cnts[i] > 0:
         stop = i
         break

   if (start == None) or (stop == None):
      return None, None, None, None, None, None
   roi  = _N.zeros(stop - start, dtype=_N.int)
   orig = _N.zeros(stop - start, dtype=_N.int)
   roi[:] = cnts[start:stop]
   
   hlfblksz = blksz/2
   blks=(stop-start)/blksz           #  full-size blocks

   hilos = _N.zeros(L, dtype=_N.int)
   hilos[0:start] = -1
   hilos[stop:]   = -1
   
   xs = range(start, start + blks*blksz, 1)
   kts= _N.zeros((2*blks)-1, _N.int)   #  full-size, but slid by hlfblksz
   
   for b in xrange((2*blks)-1):  #  slide block along by hlfblksz steps
      i = b*hlfblksz      #  i max is ((2*blks)-2)*hlfblksz
                           #  if blks == 3, blksz == 15
                           #  0 ... 4 (inclusive)

      f1 = abestfps[b, 0]
      f2 = abestfps[b, 1]
      w1 = abestfps[b, 2]
      w2 = abestfps[b, 3]

      mE = 10000    # maximum error
      bkt= -1

      #  diff of lo and hi small, or if fairly large, the weights are so
      #  skewed that it is likely that its really same hi and lo
      if (abs(f1 - f2) / (f1 + f2) < 0.25) or \
             ((abs(f1 - f2) / (f1 + f2) < 0.35) and
              ((w1 > 9*w2) or (w2 > 9*w1))):
         #  
         bkt = int(w1*f1 + w2*f2)
         if bkt == 0:   #  don't let kt be 0.  
            bkt = 1
      else:
         for kt in range(0, 3*int(f2)):
            nE = 0     #  new error
            for k in xrange(0, kt):   #  kt - 1 (inclusive)
               nE += w2 * poi_pdf(f2, k)
            for k in xrange(kt, int(f2*3)):   #  kt - 1 (inclusive)
               nE += w1 * poi_pdf(f1, k)

            if mE > nE:
               mE = nE
               bkt= kt

      kts[b] = bkt

    #  avg of kts over a range of size hlfblksz
   avgdKTs = _N.zeros(2*blks)    #  in size of hlfblksz
   avgdKTs[0] = kts[0]
   avgdKTs[2*blks-1] = kts[2*blks - 2]
   for b in xrange(1, 2*blks-1):
      avgdKTs[b] = 0.5*(kts[b - 1] + kts[b])

   #  assuming excitability slowly drifts, we linearly connect neighboring
   #  avgKTs to obtain drifting threshold
   thrshlds = _N.zeros(blks*blksz)
   thrshlds[0:hlfblksz]    = avgdKTs[0]
   thrshlds[(2*blks-1)*hlfblksz:2*blks*hlfblksz] = avgdKTs[2*blks-1]
   for b in xrange(1, 2*blks-1):
      thrshlds[b*hlfblksz:(b+1)*hlfblksz] = _N.linspace(avgdKTs[b-1], avgdKTs[b], hlfblksz)


   for b in xrange(2*blks):      #  blksz is
      for i in xrange(hlfblksz):
         hilos[start + b*hlfblksz + i] = 0
         if roi[b*hlfblksz + i] >= avgdKTs[b]:
            hilos[start + b*hlfblksz + i] = 1

   return start, start + 2 * blks, hilos, thrshlds


def rmOutlier(mdat, per=0.0, sd=0.0, col=0, topbottom=0):
   """
   Multi dimensional data.  Check column col for top and bottom per percent of data to remove
   per is between 0 and 1
   topbottom  if using percent, if topbottom == 0, remove top and bottom per percent of data.  if topbottom == -1(1), remove bottom(top) per percent only 
   stddev  how many stddevs
   """
   if len(mdat.shape) == 1:
      mdat = mdat.reshape((len(mdat), 1))
   N, D = mdat.shape

   if per != 0.0:
      scol = _N.sort(mdat[:, col])   #  sorted col
      kL   = int(per*len(scol))
      if kL == 0:
         kL = 1
      if (topbottom <= 0):
         rmLo = scol[0:kL].tolist()
      if (topbottom >= 0):
         rmHi = scol[N-kL:N].tolist()
      if topbottom < 0:
         rmHi = []
         rmdat = _N.zeros((N-kL, D), dtype=mdat.dtype)
      elif topbottom > 0:
         rmLo = []
         rmdat = _N.zeros((N-kL, D), dtype=mdat.dtype)
      else:
         rmdat = _N.zeros((N-2*kL, D), dtype=mdat.dtype)

   elif sd != 0.0:
      u    = _N.mean(mdat[:, col])
      s    = _N.std(mdat[:, col])

      rmLo = []
      rmHi = []
      for n in xrange(N):
         if mdat[n, col] > u + sd*s:
            rmHi.append(mdat[n, col])
         elif mdat[n, col] < u - sd*s:
            rmLo.append(mdat[n, col])

      rmdat = _N.zeros((N-len(rmHi)-len(rmLo), D), dtype=mdat.dtype)

   n = 0
   for i in xrange(N):
      rmFromLo = False
      rmFromHi = False
      try:
         #  don't include in rmdat
         ind = rmLo.index(mdat[i, col])
         rmLo.pop(ind)
         rmFromLo = True
      except ValueError:
         rmFromLo = False
      try:
         #  don't include in rmdat
         ind = rmHi.index(mdat[i, col])
         rmHi.pop(ind)
         rmFromHi = True
      except ValueError:
         rmFromHi = False

      if (not rmFromLo) and (not rmFromHi):
         #  didn't find it in rmLo.  Include in rmdat
         rmdat[n, :] = mdat[i, :]
         n += 1

   if D == 1:
      if per != 0.0:
         if topbottom == 0:
            rmdat = rmdat.reshape((N-2*kL, ))
         else:
            rmdat = rmdat.reshape((N-kL, ))
      else:
         rmdat = rmdat.reshape((rmdat.shape[0], ))

   return rmdat

def twostate(bestfps, pRare):
   """
   given output from 2-state EM fit (hi/lo Hzs, hi/lo weights) and a rarity
   threshold, tell me whether this set of parameters can be considered 2-state
   or not.  rarity threshold ~ 0.05 makes it so that rare outliers, which
   make EM fit send hi Hz to large values, are not taken too seriously.  
   this is my pre-conceived notion that it is an outlier value, and not a very
   rare two-state transition
   """
   segs = bestfps.shape[0]
   ans  = _N.zeros(segs, dtype=int)
   for l in xrange(segs):
      minp       = min(bestfps[l, 2], bestfps[l, 3])
      if (bestfps[l, 1] < 3):
         r      = 4
      elif (bestfps[l, 1] < 8):
         r      = 3
      else:
         r      = 2
      if (bestfps[l, 0]*r < bestfps[l, 1]) and (minp > pRare):
         ans[l] = 1
      else:
         ans[l] = 0
   return ans
