
# Leah Tomberg
import numpy

#1 Nevilles
# Declaring stuff we have
x = [3.6, 3.8, 3.9]
fx = [1.675, 1.436, 1.318]
interpoint = 3.7
def neville(x, f_x, xvalue):
    n = len(x)
    # Table for middle results
    table = [[0] * n for _ in range(n)]
    for i in range(n):
        table[i][i] = fx[i]
    # Putting vals from formula in table
    for j in range(1, n):
        for i in range(n - j):
            table[i][i + j] = ((interpoint - x[i + j]) * table[i][i + j - 1] -
                           (interpoint - x[i]) * table[i + 1][i + j]) / (x[i] -
                                                                         x[i + j])
    return table[0][n - 1]
interVal = neville(x, fx, interpoint)
print(interVal);
print("\n")


#2 Newton Forward
xi=[7.2,7.4,7.5,7.6]
fxi=[23.5492,25.3913,26.8224,27.4589]
# Degree 1 (f[x0,x1])
degree1=(fxi[1]-fxi[0])/(xi[1]-xi[0])
print(degree1)
# val under degree 1 (f[x1,x2])
underDeg1=(fxi[2]-fxi[1])/(xi[2]-xi[1])
# Degree 2
degree2=(underDeg1-degree1)/(xi[2]-xi[0])
print(degree2)
# second val under degree 1 (f[x2,x3])
secondUnderDeg1=(fxi[3]-fxi[2])/(xi[3]-xi[2])
# Val under degree 2 (f[x1,x2,x3])
underDeg2=(secondUnderDeg1-underDeg1)/(xi[3]-xi[1])
# Degree 3
degree3=(underDeg2-degree2)/(xi[3]-xi[0])
print(degree3)
print("\n")


#3 approx f(7.3)
# p1 approx
p1Approx=fxi[0]+((degree1)*(7.3-xi[0]))
# p2 approx
p2Approx=p1Approx+((degree2)*(7.3-xi[0])*(7.3-xi[1]))
#p3 approx
p3Approx=p2Approx+((degree3)*(7.3-xi[0])*(7.3-xi[1])*(7.3-xi[2]))
print(p3Approx)
print("\n")


#4 Divided difference method
xh=[3.6000000000,3.6000000000,3.8000000000,
    3.8000000000,3.9000000000,3.9000000000]
fxh=[1.6750000000,1.6750000000,1.4360000000,
     1.4360000000,1.3180000000,1.3180000000]
fprime=[-1.195,-1.188,-1.182]
# 3rd row, 3rd col
fz1z2=(fxh[2]-fxh[1])/(xh[2]-xh[1])
# 5th row, 3rd col
fz3z4=(fxh[4]-fxh[3])/(xh[4]-xh[3])
# 3rd row, 4th col
fz0z1z2=(fz1z2-fprime[0])/(xh[2]-xh[0])
# 4th row, 4th col
fz1z2z3=(fprime[1]-fz1z2)/(xh[3]-xh[1])
# 5th row, 4th col
fz2z3z4=(fz3z4-fprime[1])/(xh[4]-xh[2])
# 6th row, 4th col
fz3z4z5=(fprime[2]-fz3z4)/(xh[5]-xh[3])
# 4th row, 5th col
f0123=(fz1z2z3-fz0z1z2)/(xh[3]-xh[0])
# 5th row, 5th col
f1234=(fz2z3z4-fz1z2z3)/(xh[4]-xh[1])
# 6th row, 5th col
f2345=(fz3z4z5-fz2z3z4)/(xh[5]-xh[2])
# Formatting
numpy.set_printoptions(precision=20, suppress=True, linewidth=300)
Answer= numpy.array([[xh[0],fxh[0],0.0000000000,0.0000000000,0.0000000000],
                  [xh[1],fxh[1],fprime[0],0.0000000000,0.0000000000],
                [xh[2],fxh[2],fz1z2,fz0z1z2,0.0000000000],
                     [xh[3],fxh[3],fprime[1],fz1z2z3,f0123],
                     [xh[4],fxh[4],fz3z4,fz2z3z4,f1234],
                     [xh[5],fxh[5],fprime[2],fz3z4z5,f2345]
                  ])

print(Answer)
print("\n")

#5 Cubic spline interpolation
x5=numpy.array([2, 5, 8, 10])
f5=numpy.array([3, 5, 7, 9])
length=len(x5)
h=numpy.zeros(length)
A=numpy.zeros((length,length))
for i in range(length-1):
    h[i]=x5[i+1]-x5[i]
for i in range(1, length - 1):
    A[i, i + 1]=h[i]
# Diagonal
for i in range(1, length - 1):
    A[i,i]=2*(h[i-1]+h[i])
    A[i, i - 1]=h[i - 1]
# Making sure our bounds are all good
A[0,0]=1
A[length-1,length-1]=1
print(A)
# b stuff
b=numpy.zeros(length)
for i in range(1, length - 1):
    b[i]=3*( (f5[i + 1] - f5[i]) / h[i] - (f5[i] - f5[i - 1]) / h[i - 1] )
print(b)
xsolve=numpy.linalg.solve(A,b)
print(xsolve)
