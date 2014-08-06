import numpy as _N

class sampleARcoeff:
    pctls  = None
    cumcls = None
    k      = None

    def __init__(self, k, cls, pctls):
        self.k     = k
        self.cls   = cls
        self.pctls = pctls

        self.cumcls= _N.empty(self.k + 1)      #  cumulative cls

        self.cumcls[0]  = self.cls[0]
        for ik in xrange(1, self.k + 1):
            self.cumcls[ik] = _N.sum(self.cls[0:ik + 1])
        self.cumcls /= _N.sum(self.cls)
    
    def draw(self):
        r  = _N.random.rand()   #  for 
        ik = self.k
        nReal= 0
        while ik >= 1:
            if r >= self.cumcls[ik-1]:
                nReal = ik
                break
            ik -= 1

        nImg = self.k - nReal    #

        amps = _N.empty(nImg/2 + nReal)

        dtyp = _N.float
        if nImg > 0:
            dtyp = _N.complex
            phzs = _N.empty(nImg/2 + nReal, dtype=dtyp) * 2 * _N.pi
        else:
            phzs = _N.zeros(nReal, dtype=dtyp)

        #  For the imaginary roots
        for rr in xrange(nImg / 2):
            ###############  r
            ok = False
            while not ok:
                r = _N.random.rand()
    
                #  For the class nReal, first work with the imaginary roots, r
                end = len(self.pctls[nReal][1][0][rr][:, 1])
                ind = _N.searchsorted(self.pctls[nReal][1][0][rr][:, 1], r)
                #  For the class nReal, first work with the imaginary roots, theta
                if ind == end:
                    ind = end - 1
                amps[rr] = self.pctls[nReal][1][0][rr][ind, 0]
                if (amps[rr] <= amps[rr-1]) or (rr == 0):
                    ok = True

            ###############  theta
            r = _N.random.rand()
    
            #  For the class nReal, first work with the imaginary roots, theta
            end = len(self.pctls[nReal][1][1][rr][:, 1])
            ind = _N.searchsorted(self.pctls[nReal][1][1][rr][:, 1], r)
            #  For the class nReal, first work with the imaginary roots, theta
            if ind == end:
                ind = end - 1
            phzs[rr] = self.pctls[nReal][1][1][rr][ind, 0]

        #  For the real roots
        ofs = nImg/2 
        for rr in xrange(nImg / 2, nImg/2 + nReal):
            ###############  r
            ok = False
            while not ok:
                r = _N.random.rand()

                #  For the class nReal, first work with the imaginary roots, r
                end = len(self.pctls[nReal][0][0][rr-ofs][:, 1])
                ind = _N.searchsorted(self.pctls[nReal][0][0][rr-ofs][:, 1], r)
                #  For the class nReal, first work with the imaginary roots, theta
                if ind == end:
                    ind = end - 1
                amps[rr] = self.pctls[nReal][0][0][rr-ofs][ind, 0]
                if (amps[rr] <= amps[rr-1]) or (rr == ofs):
                    ok = True


            ###############  theta
            r = _N.random.rand()
    
            #  For the class nReal, first work with the imaginary roots, theta
            end = len(self.pctls[nReal][0][1][rr-ofs][:, 1])
            ind = _N.searchsorted(self.pctls[nReal][0][1][rr-ofs][:, 1], r)
            #  For the class nReal, first work with the imaginary roots, theta
            if ind == end:
                ind = end - 1
            phzs[rr] = self.pctls[nReal][0][1][rr-ofs][ind, 0]

        rs = _N.empty(self.k,       dtype=dtyp)
        A  = _N.zeros((self.k, self.k),  dtype=dtyp)
        B  = _N.empty((self.k, 1),  dtype=dtyp)
        F  = _N.empty((self.k, 1),  dtype=dtyp)

        for n in xrange(nImg/2):
            ph        = phzs[n]
            #  TRY SQRT
            rs[2*n]   = amps[n] * (_N.cos(ph) + 1j * _N.sin(ph))
            rs[2*n+1] = rs[2*n].real - 1j*rs[2*n].imag 
        for n in xrange(nReal):
            if phzs[nImg/2 + n] == 0:
                rs[nImg + n]   = amps[nImg/2 + n]
            else:
                rs[nImg + n]   = amps[nImg/2 + n] * -1

        for c in xrange(self.k-1):
            A[:, c] = rs**(self.k-c-1)
        A[:, self.k-1] = 1
        B[:, 0] = rs**self.k

        try:
            F = _N.linalg.solve(A, B)
        except _N.linalg.linalg.LinAlgError:
            print A, B
            F = _N.zeros((self.k, 1))
        return F.real

    def drawEG(self, nReal):
        nImg = self.k - nReal    #

        self.pctls[nReal]
        amps = _N.empty(nImg/2 + nReal)

        dtyp = _N.float
        if nImg > 0:
            dtyp = _N.complex
            phzs = _N.empty(nImg/2 + nReal, dtype=dtyp) * 2 * _N.pi
        else:
            phzs = _N.zeros(nReal, dtype=dtyp)

        #  For the imaginary roots
        for rr in xrange(nImg / 2):
            ###############  r
            r = _N.random.rand()
    
            #  For the class nReal, first work with the imaginary roots, r
            end = len(self.pctls[nReal][1][0][:, 1])
            ind = _N.searchsorted(self.pctls[nReal][1][0][:, 1], r)
            #  For the class nReal, first work with the imaginary roots, theta
            if ind == end:
                ind = end - 1
            amps[rr] = self.pctls[nReal][1][0][ind, 0]

            ###############  theta
            r = _N.random.rand()
    
            #  For the class nReal, first work with the imaginary roots, theta
            end = len(self.pctls[nReal][1][1][:, 1])
            ind = _N.searchsorted(self.pctls[nReal][1][1][:, 1], r)
            #  For the class nReal, first work with the imaginary roots, theta
            if ind == end:
                ind = end - 1
            phzs[rr] = self.pctls[nReal][1][1][ind, 0]

        #  For the real roots
        for rr in xrange(nImg / 2, nImg/2 + nReal):
            ###############  r
            r = _N.random.rand()
    
            #  For the class nReal, first work with the imaginary roots, r
            end = len(self.pctls[nReal][0][0][:, 1])
            ind = _N.searchsorted(self.pctls[nReal][0][0][:, 1], r)
            #  For the class nReal, first work with the imaginary roots, theta
            if ind == end:
                ind = end - 1
            amps[rr] = self.pctls[nReal][0][0][ind, 0]

            ###############  theta
            phzs[rr] = 0

        return amps, phzs

    def test(self, rs):  # rs Eigenvalues of F
        dtyp = rs.dtype
        A  = _N.zeros((self.k, self.k),  dtype=dtyp)
        B  = _N.empty((self.k, 1),  dtype=dtyp)
        F  = _N.empty((self.k, 1),  dtype=dtyp)

        for c in xrange(self.k-1):
            A[:, c] = rs**(self.k-c-1)
        A[:, self.k-1] = 1
        B[:, 0] = rs**self.k

        print A
        print B
        try:
            F = _N.linalg.solve(A, B)
            print "F"
            print F
        except _N.linalg.linalg.LinAlgError:
            print A, B
            F = _N.zeros((self.k, 1))
        return F.T.real
