import numpy          as _N
import kfardat        as _kfardat
import kfARlib        as _kfar
import LogitWrapper   as lw
import scipy.stats    as _ss

import warnings
warnings.filterwarnings("error")


def mcmcMixAR1(burn, NMC, y, x, zT, nStates=2, r=40, n=400, model="binomial", set=10):
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

    u   = _N.array([-2.1972245773362191, -0.40546510810816427])

    # if set==13:
    #     u   = _N.array([-1.3862943611198906, 0])
    # if set==11:
    #     u   = _N.array([-2.5866893440979424, -1.5163474893680886])
    # if set==10:
    #     u   = _N.array([-2.5866893440979424, -1.0986122886681098])

    F0  = 0.92
    q2  = 0.025
    x00 = u_x00 + _N.sqrt(s2_x00)*_N.random.rand()
    V00 = B_V00*_ss.invgamma.rvs(a_V00)

    m   = _N.random.dirichlet(alp)

    smpld_params[0, :] = (F0, q2, x00, V00)
    for i in xrange(nStates):
        smpld_params[0,i + 4] = m[i]
        smpld_params[0,i + 6] = u[i]

    smpx = x
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

        print "mean ws   %.3f" % _N.mean(ws)
        print "mean kp   %.3f" % _N.mean(kp)
        #  generate latent AR state
        _d.f_x[0, 0, 0]     = x00
        _d.f_V[0, 0, 0]     = V00
        _d.y[:]             = kp/ws - us
        _d.Rv[:] =1 / ws   #  time dependent noise

        print _N.mean(kp/ws - us)
        print _N.mean(us)
    return _d
