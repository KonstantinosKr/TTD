# simplified from http://paulbourke.net/miscellaneous/gausselim/
n = 2
print 'REAL b[%d][%d], tmp;' % (n,n)
print 'memset (b, 0, sizeof(b));'
for i in range (0,n):
  print 'b[%d][%d] = 1.0;' % (i,i)

for i in range (0,n):
  for j in range (i+1,n):
    print 'tmp = a[%d]/a[%d];' % (i*n+j,i*n+i)
    for k in range (n-1,i-1,-1):
      print 'a[%d] -= a[%d] * tmp;' % (k*n+j,k*n+i)
    for k in range (0,n):
      print 'b[%d][%d] -= b[%d][%d] * tmp;' % (k,j,k,i)

for i in range (0,n):
  for j in range (n-1,-1,-1):
    print 'tmp = 0.0;'
    for k in range (j+1, n):
      print 'tmp += a[%d] * x[%d];' % (k*n+j,i*n+k)
    print 'x[%d] = (b[%d][%d] - tmp) / a[%d];' % (i*n+j,i,j,j*n+j)
