import numpy as _N
import scipy.stats as _ss
import kstat as _ks

"""  
c functions
"""
cdef extern from "math.h":
    double exp(double)
    double sqrt(double)
    double log(double)
    double abs(double)

def gaussianEM_2comp(x, int ICs=100, double cond=0.01):
    """
    fit data cts with mixture Poisson distribution using EM
    EM depends on initial conditions.  ICs is # of initial conditions to try
    """
    cdef int N       = len(x)
    cdef int K       = 2
    pik     = _N.zeros(K)   #  weights
    uk      = _N.zeros(K)   
    sk      = _N.zeros(K)   #  std deviation
    vk      = _N.zeros(K)   #  variance
    bestu   = _N.zeros(K)
    bests   = _N.zeros(K)
    bestp   = _N.zeros(K)
    cdef double twopi = 2 * _N.pi

    s_u     = _N.mean(x)
    s_s     = _N.std(x)       
    #  need these params to check if updated component variances too small
    minDat = min(x) 
    maxDat = max(x)
    avgDist= (maxDat - minDat) / float(10*N)

    maxLL = -10000000.

    cdef int ic

    Nk      = _N.zeros(K)
    resp    = _N.zeros((N, K))

    nc      = _N.zeros(K)    # std. 
    sk2     = _N.zeros(K)    # 2*std^2

    rs      = _N.zeros(4+1)    #  place holders, random #s call at once

    #  likelihood is 1/sqrt(2 pi sig**2) * exp(-(x-u)**2/(2*sig**2))
    #  log likelihood is
    #  -0.5*log(2*pi*sig**2) + -(x-u)**2/(2*sig**2)
    llGaussian = 0
    gu      = _N.mean(x)
    sd      = _N.std(x)
    sd2     = sd*sd
    for n from 0 <= n < N:
        llGaussian += -0.5*log(twopi*sd2) - ((x[n]-gu)*(x[n]-gu))/(2*sd2)

    for ic from 0 <= ic < ICs:
        rs[:]   = _N.random.rand(5)
        pik[0]  = rs[4]
        pik[1]  = 1 - pik[0]
        uk[0] = s_u * (2*rs[0])
        uk[1] = s_u * (2*rs[1])
        sk[0] = s_s * (2*rs[2])
        sk[1] = s_s * (2*rs[3])
        vk[0] = sk[0]*sk[0]
        vk[1] = sk[1]*sk[1]

        Nk[:]         = 0.
        resp[:, :]    = 0.

        llo     = -10000000000.

        while True:
            noprob = True
            #  build responsibility array
#            rvs = []

            for k from 0 <= k < K:
                sk2[k]= 2*sk[k]*sk[k]
                nc[k] = 1/sqrt(twopi*sk2[k]*0.5)

            bCollapsed = False
            for n from 0 <= n < N:
                bot = 0

                #  bot is near 0 if component std. devs too small
                for k from 0 <= k < K:
                    bot += pik[k] * nc[k] * exp(-((x[n] - uk[k])*(x[n] - uk[k])) / sk2[k])

                if bot == 0.:
                    print "collapsed"
                    bCollapsed = True
                    resp[:, :] = 0
                if not bCollapsed:  #  otherwise get error
                    for k from 0 <= k < K:
                        resp[n, k] = (pik[k] * nc[k] * exp(-((x[n] - uk[k])*(x[n] - uk[k])) / sk2[k])) / bot
            #  build responsibility array
            if not bCollapsed:
                for k from 0 <= k < K:
                    Nk[k] = 0
                    for n from 0 <= n < N:
                        Nk[k] += resp[n, k]
                    if Nk[k] == 0.:   #  happens if sk2[k] very small
                        bCollapsed = True
            #  need new rvs
            #  update params
            if not bCollapsed:
                for k from 0 <= k < K:
                    vk[k] = 0
                    for n from 0 <= n < N:
                        vk[k] += (x[n] - uk[k])* (x[n] - uk[k])* resp[n, k]
                    vk[k] /= Nk[k]
                    sk[k] = sqrt(vk[k])
                    if sk[k] == 0.:
                        bCollapsed = True

            if not bCollapsed:
                for k from 0 <= k < K:
                    uk[k] = 0
                    for n from 0 <= n < N:
                        uk[k] += x[n] * resp[n, k]
                    uk[k] /= Nk[k]

                    pik[k] = Nk[k] / float(N)

                #  new norm parameters, since sk array updated
                for k from 0 <= k < K:
                    sk2[k]= 2*sk[k]*sk[k]
                    nc[k] = 1/sqrt(twopi*sk2[k]*0.5)

                #  log likelihood
                ll   = 0
                for n from 0 <= n < N:
                    logarg = 0
                    for k from 0 <= k < K:   #  2 components
                        logarg += pik[k] * nc[k] * exp(-((x[n] - uk[k])*(x[n] - uk[k])) / sk2[k])
                    ll += log(logarg)
                if _N.isnan(ll) or _N.isinf(ll):
                    print "nan or inf    problem"

                if (abs(llo - ll) / abs(llo + ll)) < cond:
                    break      #  stationary ll hit.  llo is log-lik of prev. iter
            else:
                #  reinitialize, hit floating point error
                ll = -1000000.
                print "hit float error, reinit params, kstat.gaussianEM_2comp()"
                noprob = False
                rs[:]   = _N.random.rand(5)
                pik[0]  = rs[4]
                pik[1]  = 1 - pik[0]
                uk[0] = s_u * (2*rs[0])
                uk[1] = s_u * (2*rs[1])
                sk[0] = s_s * (2*rs[2])
                sk[1] = s_s * (2*rs[3])
                vk[0] = sk[0]*sk[0]
                vk[1] = sk[1]*sk[1]

            if ll > maxLL:
                #  update best fit parameters
                maxLL = ll        #  component dists. wide enuff, not singularity
                bestu[0] = uk[0]
                bestu[1] = uk[1]
                bests[0] = sk[0]
                bests[1] = sk[1]
                bestp[0] = pik[0]
                bestp[1] = pik[1]
            if noprob:
                llo = ll
        
    #  I want lamk[0] to be the lower one
    if bestu[0] > bestu[1]:
        temp    = bestu[0]
        bestu[0] = bestu[1]
        bestu[1] = temp
        temp    = bestp[0]
        bestp[0]  = bestp[1]
        bestp[1]  = temp
        temp    = bests[0]
        bests[0]  = bests[1]
        bests[1]  = temp

    return bestu, bests, bestp, maxLL, llGaussian

def poissonEM_2comp(cts, int ICs=20):
    """
    fit data cts with mixture Poisson distribution using EM
    EM depends on initial conditions.  ICs is # of initial conditions to try
    """
    N       = len(cts)
    K       = 2
    pik     = _N.zeros(K)
    lamk    = _N.zeros(K)
    bestf   = _N.zeros(K)
    bestp   = _N.zeros(K)
    ICs     = 25

    maxLL = -10000000.
    uc      = _N.mean(cts)

    Nk      = _N.zeros(K)
    resp    = _N.zeros((N, K))

    logfact = _N.zeros(3*max(cts))
    rs      = _N.zeros(2+1)    #  place holders, random #s call at once

    #  log(N!) = log(N) + log(N-1) + ... + log(1)
    logfact[0] = 0
    for c from 1 <= c < 3*max(cts):
        logfact[c] = log(c) + logfact[c - 1]

    #  likelihood is lam^k exp(-lam) / fact(k)
    llPoisson = 0
    for n from 0 <= n < N:
        lpmf = cts[n]*log(uc) - uc - logfact[cts[n]]
        llPoisson += lpmf
            
    for ic from 0 <= ic < ICs:
        rs[:]   = _N.random.rand(3)
        pik[0]  = rs[2]
        pik[1]  = 1 - pik[1]
        lamk[0] = uc * 2*rs[0]
        lamk[1] = uc * 2*rs[1]

        Nk[:]         = 0.
        resp[:, :]    = 0.

        llo     = -10000000000.

        while True:
            #  build responsibility array
            for n from 0 <= n < N:
                bot = 0

                #  lam^c * exp(-lam) / c!
                for k from 0 <= k < K:
                    lpmf = (cts[n]*log(lamk[k]) - lamk[k]) - logfact[cts[n]]
                    bot += pik[k] * exp(lpmf)

                for k from 0 <= k < K:
                    lpmf = (cts[n]*log(lamk[k]) - lamk[k]) - logfact[cts[n]]
                    resp[n, k] = (pik[k] * exp(lpmf)) / bot
            #  build responsibility array
            for k from 0 <= k < K:
                Nk[k] = 0
                for n from 0 <= n < N:
                    Nk[k] += resp[n, k]

            #  update params
            for k from 0 <= k < K:
                lamk[k] = 0
                for n from 0 <= n < N:
                    lamk[k] += cts[n] * resp[n, k]
                lamk[k] /= Nk[k]

                pik[k] = Nk[k] / float(N)

            #  log likelihood
            ll   = 0
            for n from 0 <= n < N:
                logarg = 0
                for k from 0 <= k < K:
                    lpmf = cts[n]*log(lamk[k]) - lamk[k] - logfact[cts[n]]
                    logarg += pik[k] * exp(lpmf)

                ll += log(logarg)

            if (abs(llo - ll) / abs(_N.mean(llo + ll))) < 0.001:
                break

            noprob = True
            if _N.isinf(ll) or _N.isnan(ll):
                #  reinitialize, hit floating point error
#                print "hit float error, reinit params, kstat.poissonEM_2comp()"
                noprob = False
                rs[:]   = _N.random.rand(3)
                pik[0]  = rs[2]
                pik[1]  = 1 - pik[1]
                lamk[0] = uc * 2*rs[0]
                lamk[1] = uc * 2*rs[1]
            if noprob:
                llo = ll    #  last log likelihood.  each step must make smaller
            if (ll > maxLL) and noprob:
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

    return bestf, bestp, maxLL, llPoisson

def shuffle(arr):
    L = len(arr)
    t = arr[:]
    shuffled = []
    for l from 0 <= l < L:
        elem = t.pop(int(L*_N.random.rand()))
        shuffled.append(elem)
        L -= 1
    return shuffled

def shuffleThsInds(arr, inds):
    """
    instead of shuffling all elements in list, shuffle only between the indices given by inds
    """
    mtbl = []
    for l from 0 <= l < len(inds):
        mtbl.append(arr[inds[l]])  #  mutable elements

    shmtbl = shuffle(mtbl)   
    shuffled = arr[:]

    tot = 0
    for l from 0 <= l < len(inds):
        shuffled[inds[l]] = shmtbl[l]  #  mutable
    return shuffled
