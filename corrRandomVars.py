import numpy as _N

C   = 4;     # components of the noise
N   = 10000;  # data points

commonP = _N.random.normal(0, 1, N);
mu      = _N.random.rand(C);
sampMU  = _N.zeros(C);
x       = _N.zeros([C, N]);
z       = _N.zeros([C, N]);
D       = _N.random.rand(C);  # relative strength of the independent portion

for c in _N.arange(0, C):
    indep     = _N.random.normal(0, D[c], N);
    x[c] = mu[c] + commonP + indep;
    sampMU[c] = _N.mean(x[c]);

    # note that if we make D[c] too large, the mean of noises[c] will need
    # to be sampled many times before approaching mu[c].  fluctuation goes
    # as D / sqrt(N), so balance between D and N important

SIG = _N.ma.cov(x)      # SIG.data contains covariance matrix
sqSIG = db(SIG.data);   # inverse sqrt of covariance matrix

for n in _N.arange(0, N):
    # transform variables
    z[:,n] = _N.dot(_N.linalg.inv(sqSIG), x[:,n] - sampMU);

zSIG = _N.ma.cov(z)#    //  diagonal or unit variance w/ non-zero off diagonal?
