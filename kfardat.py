import numpy as _N
import os as _os

class KFARDat:
    N     = None;    #  # of data for each trial
    Ns    = None;    #  # of data for each trial
    TR    = None;    #  # of trials
    #  KF params
    p_x   = None;    p_V   = None;    p_Vi  = None
    f_x   = None;    f_V   = None
    s_x   = None;    s_V   = None
    J     = None
    Ik    = None

    #  params for AR state
    k     = None;
    ks    = None;
    F     = None
    Fs    = None
    Q     = None;
    G     = None

    #  generative and generated data
    x     = None;

    #  original spks
    dN    = None

    def __init__(self, *args):
        self.TR      = args[0]
        self.N       = args[1]
        self.k       = args[2]   #  if only using as single
        self.Ns      = _N.ones(self.TR, dtype=_N.int)*self.N
        self.ks      = _N.ones(self.TR, dtype=_N.int)*self.k

        self.initKF(self.TR, self.N, self.k)

        _N.fill_diagonal(self.F[1:, 0:self.k-1], 1)
        self.G       = _N.zeros((self.k, 1))
        self.G[0, 0] = 1
        self.Q       = _N.empty(self.TR)

    def initKF(self, TR, N, k):   #  init KF arrays
        self.F       = _N.zeros((k, k))
        self.Fs      = _N.zeros((TR, k, k))
        for tr in xrange(TR):
            _N.fill_diagonal(self.Fs[tr, 1:, 0:self.k-1], 1)
        self.Ik      = _N.identity(k)
        self.IkN     = _N.tile(self.Ik, (N+1, 1, 1))

        #  need TR
        self.p_x   = _N.empty((TR, N + 1, k, 1)) #  pr_x[:, 0]  empty, not used
        self.p_x[:, 0, 0] = 0
        self.p_V   = _N.empty((TR, N + 1, k, k))
        self.p_Vi  = _N.empty((TR, N + 1, k, k))
        self.f_x   = _N.empty((TR, N + 1, k, 1))
        self.f_V   = _N.empty((TR, N + 1, k, k))
        self.s_x   = _N.empty((TR, N + 1, k, 1))
        self.s_V   = _N.empty((TR, N + 1, k, k))
        self.J     = _N.empty((TR, N + 1, k, k))   #  K[0]  empty, not used

class KFARGauObsDat(KFARDat):
    #  Gaussian observation noise model
    y     = None
    H     = None
    Il    = None
    R     = None    #  noise
    Rv    = None    #  time-dependent noise
    K     = None

    def __init__(self, *args):
        KFARDat.__init__(self, *args)
        self.H       = _N.zeros((1, self.k))          #  row vector
        self.H[0, 0] = 1
        self.K       = _N.empty((self.TR, self.N + 1, self.k, 1))
        self.y       = _N.empty((self.TR, self.N + 1))
        self.dN      = _N.empty((self.TR, self.N + 1))
        self.Rv      = _N.empty((self.TR, self.N + 1))

    def copyData(self, *args):    #  generative data
        y     = args[0]
        if self.TR == 1:
            self.dN[0, :] = y
        else:
            self.dN[:, :] = y

    def copyParams(self, *args):   #  AR params that will get updated by EM
        F0              = args[0]
        Q               = args[1]
        self.F[0, :]    = F0[:]
        for tr in xrange(self.TR):
            self.Fs[tr, 0, :]    = F0[:]
        self.Q[:]       = Q    #  scalar
