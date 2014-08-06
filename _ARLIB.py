import numpy as _N
import sampAReig as _sae
import scipy.sparse.linalg as _ssl

def generateValidAR(k):
    """
    r    = 2
    nImg = 2   #  we want oscill.
    nReal= 0

    while r < k-2:
        if k - r > 1:
            if _N.random.rand() < 0.5:
                nImg += 2
                r += 2
            else:
                nReal += 1
                r += 1
        else:
            nReal += 1
            r += 1
    """

    r    = 0
    nImg = 0   #  we want oscill.
    nReal= 0

    while r < k:
        if k - r > 1:
            if _N.random.rand() < 0.5:
                nImg += 2
                r += 2
            else:
                nReal += 1
                r += 1
        else:
            nReal += 1
            r += 1

#    amps = 0.6*_N.random.rand(nImg/2 + nReal) + 0.4
    amps = _N.random.rand(nImg/2 + nReal)

    dtyp = _N.float
    if nImg > 0:
        dtyp = _N.complex
        phzs = _N.random.rand(nImg/2) * 2 * _N.pi

    rs = _N.empty(k,       dtype=dtyp)
    A  = _N.empty((k, k),  dtype=dtyp)
    B  = _N.empty((k, 1),  dtype=dtyp)
    F  = _N.empty((k, 1),  dtype=dtyp)

    for n in xrange(nImg/2):
        ph        = phzs[n]
        rs[2*n]   = amps[n] * (_N.cos(ph) + 1j * _N.sin(ph))
        rs[2*n+1] = rs[2*n].real - 1j*rs[2*n].imag 
    for n in xrange(nReal):
        rs[nImg + n]   = amps[nImg/2 + n]

    for c in xrange(k-1):
        A[:, c] = rs**(k-c-1)
    A[:, k-1] = 1
    B[:, 0] = rs**k
    
    F = _N.linalg.solve(A, B)
    return F.real


def cmplxRoots_old(arC):
    N   = len(arC)   # polynomial degree       a_1 B + a_2 B^2
    A   = _N.zeros((N, N))
    bBdd = True
    iBdd = 1

    A[0, :] = arC
    _N.fill_diagonal(A[1:, 0:], 1)

    vals, vecs = _N.linalg.eig(A)
    vroots = _N.empty(N)

    #  stable when roots (eigenvalues) are within unit circle
    for roots in xrange(N):
        zR = 1 / vals[roots]
        vroots[roots] = (zR * zR.conj()).real
        if vroots[roots] < 1:
            bBdd = False
            iBdd = 0

    return bBdd, iBdd, vroots, vals

def cmplxRoots(arC):
    N   = len(arC)   # polynomial degree       a_1 B + a_2 B^2
    A   = _N.zeros((N, N))
    bBdd = True
    iBdd = 1

    A[0, :] = arC
    _N.fill_diagonal(A[1:, 0:], 1)

    vals, vecs = _N.linalg.eig(A)
    rootsR     = _N.empty(N)

    #  stable when roots (eigenvalues) are within unit circle
    for roots in xrange(N):
        rootsR[roots] = _N.sqrt((vals[roots] * vals[roots].conj()).real)
        if rootsR[roots] >= 1:
            bBdd = False
            iBdd = 0

    return bBdd, iBdd, rootsR, vals

def cmplxRootsPoly(arC):
    N   = len(arC)   # polynomial degree       a_1 B + a_2 B^2
    arr = _N.empty(N+1)
    arr[0:N] = arC[::-1]
    arr[N]   = -1
    bBdd = True
    iBdd = 1
    vals = _N.polynomial.polynomial.polyroots(arr)
    rootsR     = _N.empty(N)
    #  stable when roots (eigenvalues) are within unit circle
    for roots in xrange(N):
        rootsR[roots] = _N.sqrt((vals[roots] * vals[roots].conj()).real)
        if rootsR[roots] >= 1:
            bBdd = False
            iBdd = 0
    return bBdd, iBdd, rootsR, vals

def good_prior_cov_for_AR(k):
    """
    #  for high k, relatively lower mag, higher neighbor correlation appears 
    to give reasonable AR coefficients
    """
    Mag=0.15 + 0.75/k
    cor=0.8*(1 - _N.exp(-0.3*k))
    CM = dcyCovMat(k, _N.linspace(Mag, Mag/k, k), cor)
    return CM

def diffuse_cov_for_AR(k):
    """
    #  for high k, relatively lower mag, higher neighbor correlation appears 
    to give reasonable AR coefficients
    """
    Mag=0.7
    cor=0.1
    diag = _N.linspace(Mag, Mag, k)
    diag[k-1] = 0.3
    CM = dcyCovMat(k, _N.linspace(Mag, Mag, k), cor)
    return CM

def dcyCovMat(N, diagMags, corr):
    """
    serial numbers with a decaying correlation
    N         NxN matrix
    diagMag   magnitude of the diagnals
    corr      cc of neighboring elements
    """
    covM = _N.empty((N, N))
    _N.fill_diagonal(covM, diagMags)
    for l in xrange(1, N):
        for n in xrange(N - l):
            covM[n, n+l] = (0.5*(diagMags[n]+diagMags[n+l]))*corr**l
            covM[n+l, n] = covM[n, n+l]

    return covM

def MetrSampF(N, k, smpx, q2, pr_uF, pr_cvF, burn=100, Fini=None):
    """
    sample F from full conditional given sample of latent and parameters
    N        size of latent state data
    smpx     sampled latent state
    pr_uF    prior mean for F
    pr_cvF   prior cov for F
    q2       transition noise
    From paper Wolfgang Polasek, Song Jin   "From Data to Knowledge"
    minacc   minimum number of new points accepted.  to make sure sampler
             is well-mixed.  burn alone may not be an acceptable way to
             ensure this
    """
    y   = smpx[1:N, 0].reshape(N-1, 1)    #  a column vector
    xnm = smpx[0:N-1, :]                  #  a column vector

    ipr_cvF= _N.linalg.inv(pr_cvF)        #  inverse prior cov.

    ### conditional posterior moments
    iCov   = ipr_cvF + _N.dot(xnm.T, xnm)/q2   # inv conditional posterior cov.
    Cov    = _N.linalg.inv(iCov)       # conditional posterior cov.
    M  = _N.dot(Cov, _N.dot(ipr_cvF, pr_uF) + (1/q2)*_N.dot(xnm.T, y)) # conditional posterior mean

    #   initial value of F
    if Fini == None:
        bBdd= False
        while not bBdd:
            F    = _N.random.multivariate_normal(pr_uF[:, 0], pr_cvF)  # sample from prior
            F    = F.reshape(k, 1)
            bBdd, iBdd, vr, evals = cmplxRoots(F[:, 0])
    else:
        F  = _N.empty((k, 1))
        F[:, 0] = Fini

    FM  = F - M
#    aO  = -0.5*_N.dot(FM.T, _N.dot(ipr_cvF, FM))    #  arguments to exp
    aO  = -0.5*_N.dot(FM.T, _N.dot(iCov, FM))    #  arguments to exp

    rands = _N.random.rand(burn)

    acc =0
    n = 1
    while n < burn:
        Fn  = _N.random.multivariate_normal(M[:, 0], Cov)  #  sample from cond. post.
        Fn    = Fn.reshape(k, 1)
        FnM  = Fn - M
        aC  = -0.5*_N.dot(FnM.T, _N.dot(iCov, FnM))
        bBdd, iBdd, vr, evals = cmplxRoots(Fn[:, 0])

        if bBdd:
            r   = _N.exp(aC - aO)
            if rands[n] < min(r, 1):
                acc+= 1
                F = Fn
                aO = aC
        n += 1

    return acc, F[:,0]


def MetrSampF_egv(N, k, smpx, q2, sae, burn=100, Fini=None):
    """
    sample F from full conditional given sample of latent and parameters
    N        size of latent state data
    smpx     sampled latent state
    sae      sampAReig object
    q2       transition noise
    From paper Wolfgang Polasek, Song Jin   "From Data to Knowledge"
    minacc   minimum number of new points accepted.  to make sure sampler
             is well-mixed.  burn alone may not be an acceptable way to
             ensure this
    """
    y   = smpx[1:N, 0].reshape(N-1, 1)    #  a column vector
    xnm = smpx[0:N-1, :]                  #  a column vector

    ### conditional posterior moments
    iCov   = _N.dot(xnm.T, xnm)/q2     # inv conditional posterior cov.
    Cov    = _N.linalg.inv(iCov) 
    M      = _N.dot(Cov, _N.dot(xnm.T, y))/q2 # conditional posterior mean
    # print M
    # print iCov

    #   initial value of F
    if Fini == None:
        F = sae.draw()   #  returns a column vector
    else:
        Fini = Fini.reshape((k, 1))
        F  = Fini

    FM  = F - M

    #  The Fn's being generated are not uniform in AR space
    #  This non-uniformity acts as a prior?
    aO  = -0.5*_N.dot(FM.T, _N.dot(iCov, FM))[0, 0]    #  arguments to exp

    rands = _N.random.rand(burn)

    for n in xrange(burn):
        Fn  = sae.draw()
        FnM  = Fn - M

        aC  = -0.5*_N.dot(FnM.T, _N.dot(iCov, FnM))[0, 0]
#        r   = _N.exp(aO - aC)    #  or compare aO - aC with 0
        r   = aC - aO    #  or compare aO - aC with 0
#        print "%(r).3e    %(diff).3e" % {"r" : r, "diff" : (aO-aC)}
        print "%(F)s    %(r).3e"  % {"F" : str(Fn.T), "r" : r}
        if _N.log(rands[n]) < min(r, 0):  #  
            #  switch states 
            F = Fn
            aO = aC

    return F[:,0]
