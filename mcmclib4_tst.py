import numpy          as _N
import kfardat        as _kfardat
import kfARlib        as _kfar
import LogitWrapper   as lw
import scipy.stats    as _ss


#  allow for vector observations
def mcmcMixAR1_vo(burn, NMC, y, x, zT, nStates=2, nWins=1, r=40, n=400, model="binomial"):
    """
    r   = 40    Negative binomial.  Approaches Poisson for r -> inf
    n   = 400   Binomial            
    """
    ##############################################
    ##############  Storage for sampled parameters
    ##############################################
    N       = len(y) - 1
    smpld_params = _N.empty((NMC + burn, 4 + 2*nStates))  #m1, m2, u1, u2
    z   = _N.empty((NMC+burn, N+1, nStates), dtype=_N.int)   #  augmented data

    mnCt_w1= _N.mean(y[:, 0])
    mnCt_w2= _N.mean(y[:, 1])

    #  INITIAL samples
    if model=="negative binomial":
        kp_w1   = (y[:, 0] - r[0]) *0.5
        kp_w2   = (y[:, 1] - r[1]) *0.5
        p0_w1   = mnCt_w1 / (mnCt_w1 + r[0])       #  matches 1 - p of genearted
        p0_w2   = mnCt_w2 / (mnCt_w2 + r[0])       #  matches 1 - p of genearted
        rn   = r    #  length nWins
    else:
        kp_w1  = y[:, 0] - n[0]*0.5
        kp_w2  = y[:, 1] - n[1]*0.5
        p0_w1  = mnCt_w1 / float(n[0])       #  matches 1 - p of genearted
        p0_w2  = mnCt_w2 / float(n[1])       #  matches 1 - p of genearted
        rn   = n    #  length nWins
    u0_w1  = _N.log(p0_w1 / (1 - p0_w1))    #  -1*u generated
    u0_w2  = _N.log(p0_w2 / (1 - p0_w2))    #  -1*u generated

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
    u_u        = _N.empty(nStates + nWins)
    s2_u       = _N.zeros((nStates + nWins, nStates + nWins))
    # (win1 s1)    (win1 s2)    (win2 s1)    (win2 s2)
    u_u[:]     = (u0_w1*1.2, u0_w1*0.8, u0_w2*1.2, u0_w2*0.8)
    _N.fill_diagonal(s2_u, [0.5, 0.5, 0.5, 0.5])
    #  q2  --  Inverse Gamma prior
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
    _d.copyData(y[:, 0], y[:, 0])   #  dummy data copied

    u_w1   = _N.array([-2.5866893440979424, -1.0986122886681098])
    u_w2   = _N.array([-2.5866893440979424, -1.0986122886681098])

#    u_w1   = _N.random.multivariate_normal(u_u[0:nStates], s2_u[0:nStates, 0:nStates])
#    u_w2   = _N.random.multivariate_normal(u_u[nStates:2*nStates], s2_u[nStates:nStates*2, nStates:nStates*2])

    F0  = ((b_F0) - (a_F0)) * _N.random.rand() + a_F0
    F0 = 0.92
    q2  = B_q2*_ss.invgamma.rvs(a_q2)
    q2  = 0.015
    x00 = u_x00 + _N.sqrt(s2_x00)*_N.random.rand()
    V00 = B_V00*_ss.invgamma.rvs(a_V00)

#    m   = _N.random.dirichlet(alp)
    m   = _N.zeros(nStates)
    m[0] = _N.sum(1-zT) / float(N+1)
    m[1] = 1 - m[0]
    smp_F        = _N.zeros(NMC + burn)
    smp_q2       = _N.zeros(NMC + burn)
    smp_u        = _N.zeros((NMC + burn, nWins, nStates))   #  uL_w1, uH_w1, uL_w2, uH_w2, ....
    smp_m        = _N.zeros((NMC + burn, nStates))

    smpx = _N.zeros(N + 1)   #  start at 0 + u
    smpx[:] = x[:]   #  start at 0 + u
    Bsmpx= _N.zeros((NMC, N + 1))

    ws_w1 = lw.rpg_devroye(rn[0], smpx + u0_w1, num=(N + 1))
    ws_w2 = lw.rpg_devroye(rn[1], smpx + u0_w2, num=(N + 1))

    trm_w1  = _N.empty(nStates)
    trm_w2  = _N.empty(nStates)

    for it in xrange(1, NMC+burn):
        if (it % 50) == 0:
            print it
        #  generate latent zs.  Depends on Xs and PG latents
        kw_w1  = kp_w1 / ws_w1
        kw_w2  = kp_w2 / ws_w2
        rnds =_N.random.rand(N+1)
        z[it, :, 0] = 1-zT
        z[it, :, 1] = zT

        #  generate PG latents.  Depends on Xs and us, zs.  us1 us2 
        us_w1 = _N.dot(z[it, :, :], u_w1)   #  either low or high uf
        us_w2 = _N.dot(z[it, :, :], u_w2)
        #  us_w1 is like -->  [u1  u2  u2  u2  u1  u1  u1 ...]
        ws_w1 = lw.rpg_devroye(rn[0], smpx + us_w1, num=(N + 1))
        ws_w2 = lw.rpg_devroye(rn[1], smpx + us_w2, num=(N + 1))
        _d.copyParams(_N.array([F0]), q2, _N.array([1]), 1)

        #  generate latent AR state
        _d.f_x[0, 0, 0]     = x00
        _d.f_V[0, 0, 0]     = V00
        btm      = 1 / ws_w1 + 1 / ws_w2   #  shape N x 1

        top = (kw_w1 - us_w1) / ws_w2 + (kw_w2 - us_w2) / ws_w1

        _d.y[:] = top/btm
        _d.Rv[:] =1 / (ws_w1 + ws_w2)   #  time dependent noise
        smpx = _kfar.armdl_FFBS_1itr(_d, samples=1)

        #  p3 --  samp u here

        dirArgs = _N.empty(nStates)
        for i in xrange(nStates):
            dirArgs[i] = alp[i] + _N.sum(z[it, :, i])
        m[:] = _N.random.dirichlet(dirArgs)

        # # sample u
        for st in xrange(nStates):
            #  win1 for this state
            iw = st + 0 * nStates

            A = 0.5*(1/s2_u[iw,iw] + _N.dot(ws_w1, z[it, :, st]))
            B = u_u[iw]/s2_u[iw,iw] + _N.dot(kp_w1 - ws_w1*smpx, z[it, :, st])
            u_w1[st] = B/(2*A) + _N.sqrt(1/(2*A))*_N.random.randn()
#            print "mean   u_w1[%(st)d] = %(u).3f" % {"st" : st, "u" : u_w1[st]}
            #  win2 for this state
            iw = st + 1 * nStates

            A = 0.5*(1/s2_u[iw,iw] + _N.dot(ws_w2, z[it, :, st]))
            B = u_u[iw]/s2_u[iw,iw] + _N.dot(kp_w2 - ws_w2*smpx, z[it, :, st])
            u_w2[st] = B/(2*A) + _N.sqrt(1/(2*A))*_N.random.randn()
#            print "mean   u_w2[%(st)d] = %(u).3f" % {"st" : st, "u" : u_w2[st]}
            

        # sample F0
        F0AA = _N.dot(smpx[0:-1], smpx[0:-1])
        F0BB = _N.dot(smpx[0:-1], smpx[1:])

        F0std= _N.sqrt(q2/F0AA)
        F0a, F0b  = (a_F0 - F0BB/F0AA) / F0std, (b_F0 - F0BB/F0AA) / F0std
        F0=F0BB/F0AA+F0std*_ss.truncnorm.rvs(F0a, F0b)

#        print "%(F0).4f    %(q2).4f" % {"F0" : F0, "q2" : q2}

        #####################    sample q2
        a = a_q2 + 0.5*(N+1)  #  N + 1 - 1
        rsd_stp = smpx[1:] - F0*smpx[0:-1]
        BB = B_q2 + 0.5 * _N.dot(rsd_stp, rsd_stp)
        q2 = _ss.invgamma.rvs(a, scale=BB)
#        print q2
        # #####################    sample x00
        mn  = (u_x00*V00 + s2_x00*x00) / (V00 + s2_x00)
        vr = (V00*s2_x00) / (V00 + s2_x00)
#        x00 = mn + _N.sqrt(vr)*_N.random.randn()
        #####################    sample V00
        aa = a_V00 + 0.5
        BB = B_V00 + 0.5*(smpx[0] - x00)*(smpx[0] - x00)
#        V00 = _ss.invgamma.rvs(aa, scale=BB)

        smp_F[it]       = F0
        smp_q2[it]      = q2
        smp_u[it, 0, :] = u_w1
        smp_u[it, 1, :] = u_w2
        smp_m[it, :]    = m

        if it >= burn:
            Bsmpx[it-burn, :] = smpx

    return Bsmpx, smp_F, smp_q2, smp_u, smp_m, z, _d
