import numpy          as _N
import kfardat        as _kfardat
import kfARlib        as _kfar
import LogitWrapper   as lw
import scipy.stats    as _ss

import warnings
warnings.filterwarnings("error")


def mcmcMixAR1(burn, NMC, y, nStates=2, r=40, n=400, model="binomial"):
    """
    r   = 40    Negative binomial.  Approaches Poisson for r -> inf
    n   = 400   Binomial            
    """
    ##############################################
    ##############  Storage for sampled parameters
    ##############################################
    N       = len(y) - 1
    nStates = 2
    smpld_params = _N.empty((NMC + burn, 4 + 2*nStates))  #m1, m2, u1, u2
    z   = _N.empty((NMC+burn, N+1, nStates))   #  augmented data

    mnCt= _N.mean(y)
    #  INITIAL samples
    if model=="negative binomial":
        kp   = (y - r) *0.5
        p0   = mnCt / (mnCt + r)       #  matches 1 - p of genearted
        rn   = r
    else:
        kp  = y - n*0.5
        p0   = mnCt / float(n)       #  matches 1 - p of genearted
        rn   = n
    u0  = _N.log(p0 / (1 - p0))    #  -1*u generated

    #######  PRIOR parameters
    #  F0  --  flat prior
    #a_F0         = -1
    #  I think a prior assumption of relatively narrow and high F0 range
    #  is warranted.  Small F0 is close to white noise, and as such can
    #  be confused with the independent trial-to-trial count noise.Force
    #  it to search for longer timescale correlations by setting F0 to be
    #  fairly large.
    a_F0         = -0.1             #  prior assumption: slow fluctuation
    b_F0         =  1
    #  u   --  Gaussian prior
    u_u        = _N.empty(nStates)
    s2_u       = _N.zeros((nStates, nStates))
    u_u[:]     = (u0*1.2, u0*0.8)
    _N.fill_diagonal(s2_u, [0.5, 0.5])
    #  q2  --  Inverse Gamma prior
#    pr_mn_q2     = 0.05
    pr_mn_q2     = 0.05
    a_q2         = 2
    B_q2         = (a_q2 + 1) * pr_mn_q2
    #  x0  --  Gaussian prior
    u_x00        = 0
    s2_x00       = 0.5
    #  V00 --  Inverse Gamma prior
    pr_mn_V00    = 1    #  mode of prior distribution of variance ~ 1
    a_V00        = 2
    B_V00        = (a_V00 + 1)*pr_mn_V00
    #  m1, m2 --  Dirichlet prior
    alp          = _N.ones(nStates)

    # #generate initial values of parameters
    #generate initial time series
    _d = _kfardat.KFARGauObsDat(N, 1)
#    _d.copyData(y, x_st_cnts[:, 0])
    _d.copyData(y, y)   #  dummy data copied

    u   = _N.random.multivariate_normal(u_u, s2_u)
    u   = _N.array([-2.5866893440979424, -1.0986122886681098])
#    F0  = ((b_F0) - (a_F0)) * _N.random.rand() + a_F0
    F0  = 0.92
#    q2  = B_q2*_ss.invgamma.rvs(a_q2)
    q2  = 0.015
    x00 = u_x00 + _N.sqrt(s2_x00)*_N.random.rand()
    V00 = B_V00*_ss.invgamma.rvs(a_V00)

#    m   = _N.random.dirichlet(alp)
    m   = _N.zeros(nStates)
    m[0] = zT[0] / float(N+1)
    m[1] = zT[1] / float(N+1)

    smpld_params[0, :] = (F0, q2, x00, V00)
    for i in xrange(nStates):
        smpld_params[0,i + 4] = m[i]
        smpld_params[0,i + 6] = u[i]

    smpx = _N.zeros(N + 1)   #  start at 0 + u
    Bsmpx= _N.zeros((NMC, N + 1))
    trm  = _N.empty(nStates)

    ws = lw.rpg_devroye(rn, smpx + u0, num=(N + 1))

    for it in xrange(1, NMC+burn):
        if (it % 50) == 0:
            print it
        #  generate latent zs.  Depends on Xs and PG latents
        kw  = kp / ws
        rnds =_N.random.rand(N+1)
        z[it, :, 0] = 1-zT
        z[it, :, 1] = zT


        #  generate PG latents.  Depends on Xs and us, zs
        us   = _N.dot(z[it, :, :], u)
        ws = lw.rpg_devroye(rn, smpx + us, num=(N + 1))
        _d.copyParams(_N.array([F0]), q2, _N.array([1]), 1)

        #  generate latent AR state
        _d.f_x[0, 0, 0]     = x00
        _d.f_V[0, 0, 0]     = V00
        _d.y[:]             = kp/ws - us
        _d.Rv[:] =1 / ws   #  time dependent noise
        smpx = _kfar.armdl_FFBS_1itr(_d, samples=1)

        #  p3 --  samp u here

        dirArgs = _N.empty(nStates)
        for i in xrange(nStates):
            dirArgs[i] = alp[i] + _N.sum(z[it, :, i])
        m[:] = _N.random.dirichlet(dirArgs)

        # # sample u
        for i in xrange(nStates):
            A = 0.5*(1/s2_u[i,i] + _N.dot(ws, z[it, :, i]))
            B = u_u[i]/s2_u[i,i] + _N.dot(kp - ws*smpx, z[it, :, i])
            u[i] = B/(2*A) + _N.sqrt(1/(2*A))*_N.random.randn()
            print "mean   u[%(st)d] = %(u).3f" % {"st" : i, "u" : u[i]}

        # sample F0
        F0AA = _N.dot(smpx[0:-1], smpx[0:-1])
        F0BB = _N.dot(smpx[0:-1], smpx[1:])

        F0std= _N.sqrt(q2/F0AA)
        F0a, F0b  = (a_F0 - F0BB/F0AA) / F0std, (b_F0 - F0BB/F0AA) / F0std
        F0=F0BB/F0AA+F0std*_ss.truncnorm.rvs(F0a, F0b)

        #####################    sample q2
        a = a_q2 + 0.5*(N+1)  #  N + 1 - 1
        rsd_stp = smpx[1:] - F0*smpx[0:-1]
        BB = B_q2 + 0.5 * _N.dot(rsd_stp, rsd_stp)
#        print BB / (a-1)
        q2 = _ss.invgamma.rvs(a, scale=BB)
        #####################    sample x00
        mn  = (u_x00*V00 + s2_x00*x00) / (V00 + s2_x00)
        vr = (V00*s2_x00) / (V00 + s2_x00)
        x00 = mn + _N.sqrt(vr)*_N.random.randn()
        #####################    sample V00
        aa = a_V00 + 0.5
        BB = B_V00 + 0.5*(smpx[0] - x00)*(smpx[0] - x00)
        V00 = _ss.invgamma.rvs(aa, scale=BB)

        smpld_params[it, 0:4] = (F0, q2, x00, V00)
        for i in xrange(nStates):
             smpld_params[it,i + 4] = m[i]
             smpld_params[it,i + 6] = u[i]
        if it >= burn:
            Bsmpx[it-burn, :] = smpx

    return Bsmpx, smpld_params, z

