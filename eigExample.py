A = _N.array([[1.0, 3, -1], [2, 0, 3], [1, 1, 9]])

w,v=_N.linalg.eig(A)
# column vector v[:,i] corresponds to eigenvalue w[i]

w[0] * v[:, 0] - _N.dot(A, v[:,0])   # should be nearly 0
w[1] * v[:, 1] - _N.dot(A, v[:,1])   # should be nearly 0
w[2] * v[:, 2] - _N.dot(A, v[:,2])   # should be nearly 0



S = _N.array([[1.0, 3, -1], [3, 0, 2], [-1, 2, 9]])

sw,sv=_N.linalg.eig(S)
# column vector v[:,i] corresponds to eigenvalue w[i]

sw[0] * sv[:, 0] - _N.dot(S, sv[:,0])   # should be nearly 0
sw[1] * sv[:, 1] - _N.dot(S, sv[:,1])   # should be nearly 0
sw[2] * sv[:, 2] - _N.dot(S, sv[:,2])   # should be nearly 0
