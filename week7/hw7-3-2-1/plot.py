#!/usr/loca/bin/python
from matplotlib import rc
from matplotlib import rc
import matplotlib.pyplot as plt
from numpy import *
from pylab import *

data = loadtxt("u.txt");

m = max(data[:,0])+1;
n = max(data[:,1])+1

m = int(m)
n = int(n)

X = zeros([m, n]);
Y = zeros([m, n]);
field = zeros([m, n]);

for q in arange(0,m*n,1):
        i = int(data[q,0])
        j = int(data[q,1])
        field[i,j] = data[q,4]
        x = data[q,2]
        y = data[q,3]
        X[i, j] = x
        Y[i, j] = y

plt.figure(1)
plt.contourf(X,Y,field,20,cmap=cm.jet)
plt.colorbar()
plt.title("u field",fontsize=20)

plt.show()
