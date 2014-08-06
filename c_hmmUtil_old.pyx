#  hand me transition matrix
import numpy as _N

"""  
c functions
"""
cdef extern from "math.h":
    double log(double)
    double exp(double)
    double floor(double)

def sumOfLogs(y, N):
    Os = _N.empty(N)
    Ms = _N.empty(N)
    #  sum_i p_i    each p_i small, this sum == 0 if done naively
    #  store each p_i as a log

    LL = 0
    for i from 0 <= i < N:
        Os[i] = floor(y[i])
        Ms[i] = exp(y[i] - Os[i])

    maxO = max(Os)
    M    = 0
    for i from 0 <= i < N:
        M += Ms[i] * exp((Os[i] - maxO))

    LL = log(M) + maxO
    return LL

def sumOfLogsCond(y, N, O, k):
    Os = _N.empty(N)
    Ms = _N.empty(N)
    #  sum_i p_i    each p_i small, this sum == 0 if done naively
    #  store each p_i as a log

    LL = 0
    maxO = -1000000000   #  
    for i from 0 <= i < N:
        if O[i] == k:
            Os[i] = floor(y[i])
            Ms[i] = exp(y[i] - Os[i])
            if Os[i] > maxO:
                maxO = Os[i]

    M    = 0
    for i from 0 <= i < N:
        if O[i] == k:
            M += Ms[i] * exp((Os[i] - maxO))

    LL = log(M) + maxO
    return LL

def Lfrwd(O, N, T, a, b, pi):
    cdef int i, j, t
    cdef double LL
    
    Os = _N.empty(N)
    Ms = _N.empty(N)
    al = _N.zeros((N, T))   #  al is a log probability

    #  Initialize al
    for i from 0 <= i < N:
        al[i, 0] = log(pi[i]) + log(b[i, O[0]])

    #  Calculate a new al for (t + 1)th time step 
    for t from 0 <= t < T - 1:
        for j from 0 <= j < N:
            for i from 0 <= i < N:
                L = al[i, t] + log(a[i, j])
                Os[i] = floor(L)
                Ms[i] = exp(L - Os[i])

            maxO = max(Os)
            M    = 0
            for i from 0 <= i < N:
                M += Ms[i] * exp((Os[i] - maxO))

            Lsum = log(M) + maxO

            #  log probability complete   #  sum of PRODUCT terms
            al[j, t + 1] = log(b[j, O[t + 1]]) + Lsum

    LL = sumOfLogs(al[:, T - 1], N)

    return LL, al

#  hand me transition matrix
def Lbkwd(O, N, T, a, b):
    Os = _N.empty(N)
    Ms = _N.empty(N)

    bt = _N.zeros((N, T))   # beta

    for i from 0 <= i < N:
        bt[i, T - 1] = 0

    #  Calculate a new al for t-th time step
    for t from T - 2 >= t > -1:
        for i from 0 <= i < N:
            for j from 0 <= j < N:
                L = log(a[i, j]) + log(b[i, O[t + 1]]) + bt[j, t + 1]
                Os[j] = floor(L)       #  Order of magnitude
                Ms[j] = exp(L - Os[j]) #  left over

            maxO = max(Os)
            M    = 0
            for j from 0 <= j < N:
                M += Ms[j] * exp((Os[j] - maxO))  #  summing the "left overs"

            Lsum = log(M) + maxO  #  

            #  log probability complete
            bt[i, t] = Lsum
    return bt

#  Baum-Welch
def BW(O, N):
    #  for calculating log of sums
    Os = _N.empty(N)
    Ms = _N.empty(N)

    E_pi = _N.empty(N)
    E_a  = _N.empty((N, N))
    E_b  = _N.empty((N, N))

    a  = _N.zeros((N, N))
    b  = _N.zeros((N, N))
    pi = _N.zeros(N)
    fTUD = 0.4                 # for transitions from up to down
    fTDU = 0.4                 # for transitions from up to down
    dt = 0.5 * 0.001
    fD   = 2.5
    fU   = 15.2
    a[0, 0] = 1 - fTDU * dt    # keep down state
    a[0, 1] = fTDU * dt        # down to up transition
    a[1, 0] = fTUD * dt        # up to down transition
    a[1, 1] = 1 - fTUD * dt    # keep up state
    b[0, 0] = 1 - fD * dt      # in down state, no spike
    b[0, 1] = fD * dt          # in down state, spike
    b[1, 0] = 1 - fU * dt      # in up state, no spike
    b[1, 1] = fU * dt          # in up state, spike
    pi[0]   = 0.5              # down state prior probability 
    pi[1]   = 0.5              # down state prior probability 

    T       = len(O)

    LL, al = Lfrwd(O, N, T, a, b, pi)
    bt     = Lbkwd(O, N, T, a, b)

    print "fD   %(fd).1f  fU   %(fu).1f" % {"fd" : fD, "fu" : fU}
    print "fTUD %(fd).1f  fTDU %(fu).1f" % {"fd" : fTUD, "fu" : fTDU}
    print "%.3e" % LL
    
    xi    = _N.empty((N, N, T - 1))   #  log probability

    for t from 0 <= t < T - 1:
        for i from 0 <= i < N:
            for j from 0 <= j < N:
                xi[i, j, t] = al[i, t] + log(a[i, j]) + log(b[j, O[t + 1]]) + bt[j, t + 1] - LL

    # for each time, we want gamma_t(i)
    gam = _N.empty((N, T - 1))
    occ = _N.empty(N)   #  expected occupancy of state i
    #  we have a list of log probabilities.  want log of sum of probs.
    for t from 0 <= t < T - 1:
        for i from 0 <= i < N:
            gam[i, t] = sumOfLogs(xi[i, :, t], N)

    for i from 0 <= i < N:
        occ[i] = sumOfLogs(gam[i, :], T - 1)

    for i from 0 <= i < N:
        E_pi[i] = exp(gam[i, 0])
    for i from 0 <= i < N:
        for j from 0 <= j < N:        
            E_a[i, j] = sumOfLogs(xi[i, j, :], T - 1) - occ[i]
    for i from 0 <= i < N:
        for j from 0 <= j < N:
            for k from 0 <= k < N:
                E_b[i, k] = sumOfLogsCond(gam[i, :], T - 1, O, k) - occ[i]
        
    return xi, occ, E_pi, E_a, E_b
#     print E_pi
#     return E_pi
    
