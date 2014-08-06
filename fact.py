nn = 1.

for n in xrange(1, 60):
    #  start with 1!
    nn *= n
    print "%(1).5e %(2).5e" % {"1" : _N.log(nn), "2" : (n*_N.log(n) - n)}
