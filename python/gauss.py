# simplified from http://paulbourke.net/miscellaneous/gausselim/
n = 2
print 'REAL tmp;' % (n,n)
for i in range (0,n):
  for j in range (i+1,n):
    print 'tmp = a[%d]/a[%d];' % (i*n+j,i*n+i)
    for k in range (n-1,i-1,-1):
      print 'a[%d] -= a[%d] * tmp;' % (k*n+j,k*n+i)
    print 'b[%d] -= b[%d] * tmp;' % (j,i)

for j in range (n-1,-1,-1):
  print 'tmp = 0.0;'
  for k in range (j+1, n):
    print 'tmp += a[%d] * x[%d];' % (k*n+j,k)
  print 'x[%d] = (b[%d] - tmp) / a[%d];' % (j,j,j*n+j)
