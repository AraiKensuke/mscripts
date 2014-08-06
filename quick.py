nCS = 2
nCN = 10

for n in xrange(nCS + nCN):
    if n >= nCS:
        th = (0.15 + 0.8 / (nCN - (n - nCS)))
        r  = (1 + _N.random.rand()) / (2 * _N.sqrt(nCN))
        print "%(th).2f   %(r).2f" % {"th" : th, "r" : r}

