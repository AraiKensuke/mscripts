#  given covariance matrix and list of indices, tell me what the 

def marginalPDF(cov, lndx):
    """
    cov    covariance matrix 
    lndx   list of indices to marginalize over
    """

def practice():    
    #  remove rows and columns from matrix

    A = _N.array([[ 1, -1, -2, -3, -4, -5],
                  [-1,  2,  0,  0,  0,  0],
                  [-2,  0,  3,  0,  0,  0],
                  [-3,  0,  0,  4,  0,  0],
                  [-4,  0,  0,  0,  5,  0],
                  [-5,  0,  0,  0,  0,  6]])

    lndx = _N.array([2, 4])
    ulndx= _N.array([0, 1, 3, 5])

    B = _N.empty((6, 6))
    
    B[0:2, 0:2] = A[_N.ix_(lndx,  lndx)]
    B[2:,  2:]  = A[_N.ix_(ulndx, ulndx)]
    B[0:2, 2:]  = A[_N.ix_(lndx, ulndx)]
    B[2:,  0:2] = A[_N.ix_(ulndx, lndx)]


