import time as _tm

k = 10
I = _N.identity(k)
N = 20000

f_x = _N.random.multivariate_normal(_N.zeros(k), I, size=N)
p_x = _N.zeros((N, k))

F = _N.zeros((k, k))
_N.fill_diagonal(F[1:, 0:k-1], 1)
F[0, :] = _N.random.randn(k)

t1 = _tm.time()
for n in xrange(N):
    p_x[n] = _N.dot(F, f_x[n])
t2 = _tm.time()
print p_x

F1 = _N.array(F[0, :])
t3 = _tm.time()
tmp = _N.empty(k)
for n in xrange(N):
    tmp[0]  = _N.dot(F1, f_x[n])
    tmp[1:] = f_x[n, 0:k-1]
    p_x[n] = tmp[:]
t4 = _tm.time()

# t1 = _tm.time()
# for n in xrange(N):
#     _N.dot(F, f_x[n])
# t2 = _tm.time()

print (t2-t1)
print (t4-t3)
