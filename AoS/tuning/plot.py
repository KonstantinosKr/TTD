import numpy as np
import pandas as pd
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.mlab import griddata
from scipy import stats
from mpl_toolkits.mplot3d import axes3d

pd.set_option('display.width', 1000)
pd.set_option('display.max_rows', 100)

def plotTuneFail(mthd):
    neps = int(mthd.ix[0, 7])
    nr = int(mthd.ix[0, 8])
    nlgth = int(mthd.ix[0, 9])

    xx = np.zeros((nr, neps))
    yy = np.zeros((nr, neps))
    zz = np.zeros((nr, neps))

    fig = plt.figure()
    for k in range (0,nlgth):
        m = nr*neps*k
        for i in range (0,nr):
            for j in range (0,neps):
                xx[i,j] = mthd.ix[m+i*neps+j,2]#r
                yy[i,j] = mthd.ix[m+i*neps+j,3]#eps
                zz[i,j] = mthd.ix[m+i*neps+j,4]#failre

        xx = np.log10(xx)
        yy = np.log10(yy)

        ax = fig.add_subplot(np.ceil(float(nlgth/2.)),2,k+1, projection='3d')
        surf = ax.plot_surface(xx, yy, zz, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
        #plt.title('Tuning (r x eps x failure%)');
        plt.xlabel('r parameter', fontsize=16)
        plt.ylabel('eps parameter', fontsize=16)
    plt.show()
    return;

def plotTuneIt(mthd):
    neps = int(mthd.ix[0, 7])
    nr = int(mthd.ix[0, 8])
    nlgth = int(mthd.ix[0, 9])

    xx = np.zeros((nr, neps))
    yy = np.zeros((nr, neps))
    zz = np.zeros((nr, neps))

    fig = plt.figure()
    for k in range (0,nlgth):
        m = nr*neps*k
        for i in range (0,nr):
            for j in range (0,neps):
                xx[i,j] = mthd.ix[m+i*neps+j,2]#r
                yy[i,j] = mthd.ix[m+i*neps+j,3]#eps
                zz[i,j] = mthd.ix[m+i*neps+j,5]#iterations

        xx = np.log10(xx)
        yy = np.log10(yy)

        ax = fig.add_subplot(np.ceil(float(nlgth/2.)),2,k+1, projection='3d')
        surf = ax.plot_surface(xx, yy, zz, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
        #plt.title('Tuning (r x eps x failure%)');
        plt.xlabel('r parameter', fontsize=25)
        plt.ylabel('eps parameter', fontsize=25)
    plt.show()
    return;

def plotError(mthd):
    neps = int(mthd.ix[0, 7])
    nr = int(mthd.ix[0, 8])
    nlgth = int(mthd.ix[0, 9])

    xx = np.zeros((nr, neps))
    yy = np.zeros((nr, neps))
    zz = np.zeros((nr, neps))

    fig = plt.figure()
    for k in range (0,nlgth):
        m = nr*neps*k
        for i in range (0,nr):
            for j in range (0,neps):
                xx[i,j] = mthd.ix[m+i*neps+j,2]#r
                yy[i,j] = mthd.ix[m+i*neps+j,3]#eps
                zz[i,j] = mthd.ix[m+i*neps+j,6]#error

        xx = np.log10(xx)
        yy = np.log10(yy)

        ax = fig.add_subplot(np.ceil(float(nlgth/2.)),2,k+1, projection='3d')
        surf = ax.plot_surface(xx, yy, zz, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
        plt.xlabel('r parameter', fontsize=25)
        plt.ylabel('eps parameter', fontsize=25)
    plt.show()
    return;

penalty = pd.read_csv("ptune.data")

#plotTuneFail(penalty);
plotError(penalty);
plotTuneIt(penalty);



