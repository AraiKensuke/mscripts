import scipy as _S

__gammacdf = []
__gammacdfkappa = []

def _gammaPdist(k, x):
    return k*(k*x)**(k - 1)*_N.exp(-k*x)/_S.special.gamma(k)

#  sample from Gamma dist where theta = 1/k
def _sampleUnitGammaDist(k, samples=1):
    iUseCDF = -1
    try:
        iUseCDF = __gammacdfkappa.index(k)
    except:
        print "making new CDF"
        __gammacdfkappa.append(k)
        iUseCDF = len(__gammacdfkappa)

        cv = 0
        x  = 0.
        dx = 0.001
        dx2= dx*0.5

        while cv < 0.999:
            a = dx * _gammaPdist(k, x + dx2)
            cv += a
            __gammacdf.append(a)
            
    for s in range(samples):
        
