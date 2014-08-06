import numpy as _N
import os as _os

class KFARDatMT:
    N     = None;    #  # of data for each trial
    Tr    = None;    #  # of trials
    #  KF params
    p_x   = None;    p_V   = None;    p_Vi  = None
    f_x   = None;    f_V   = None
    s_x   = None;    s_V   = None
    J     = None
    Ik    = None

    #  params for AR state
    k     = None;
    F     = None
    Q     = None;
    G     = None

    #  generative and generated data
    x     = None;

    #  Pi   (permutation matrix, Fruhwirth-Schnetter)
    Pi    = None;
    Pi0   = None;   #  Pi0   r-d x r
    Pist  = None;   #  Pi*   d   x r

    def __init__(self, *args):
        self.N            = args[0]

        k            = args[1]
        self.multiAR       = False
        self.k       = k   #  if only using as single
        self.initKF(self.N, k)
        self.Pi      = _N.identity(k)
        self.Pi0     = self.Pi[1:, :]
        self.Pist    = self.Pi[0, :]

        for i in xrange(1, k):
            self.F[i, i-1] = 1
        self.G       = _N.zeros((k, 1))
        self.G[0, 0] = 1

    def initKF(self, N, k):   #  init KF arrays
        self.F       = _N.zeros((k, k))
        self.p_x   = _N.empty((N + 1, k, 1))   #  pr_x[:, 0]  empty, not used
        self.p_x[:, 0, 0] = 0
        self.p_V   = _N.empty((N + 1, k, k))
        self.p_Vi  = _N.empty((N + 1, k, k))
        self.f_x   = _N.empty((N + 1, k, 1))
        self.f_V   = _N.empty((N + 1, k, k))
        self.s_x   = _N.empty((N + 1, k, 1))   #  pr_x[:, 0]  empty, not used
        self.s_V   = _N.empty((N + 1, k, k))
        self.J     = _N.empty((N + 1, k, k))   #  K[0]  empty, not used
        self.Ik      = _N.identity(k)
        self.IkN     = _N.tile(self.Ik, (N+1, 1, 1))

class KFARGauObsDatMT(KFARDat):
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
        self.K       = _N.empty((self.N + 1, self.k, 1))


    def copyData(self, *args):    #  generative data
        x     = args[0]
        y     = args[1]
        self.x   = _N.empty(self.N + 1)
        self.x[:]  = x[:]
        self.y   = _N.empty(self.N + 1)
        self.y[:]  = y[:]

    def copyParams(self, *args):   #  AR params that will get updated by EM
        F0              = args[0]
        Q               = args[1]
        H               = args[2]
        R               = args[3]
        self.F[0, :]    = F0[:]
        self.Q          = Q    #  scalar
        self.H[:]       = H[:]
        self.R          = R
        self.Rv         = _N.ones(self.N + 1)*R

