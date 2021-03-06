from ctypes import *
from random import shuffle, sample
import time
"""
libfunctions = cdll.LoadLibrary("./libfunctions.so")
numbers = sample(range(1000000), 99999)
shuffle(numbers)
c_numbers = (c_int * len(numbers))(*numbers)


start = time.time()
libfunctions.merge_sort(byref(c_numbers), len(numbers))
finish = time.time()
print("C: " + str(finish - start))


"""
def mainW(cost):
    def costNoSt(timePoint):
        c1 = np.asarray(1./100000.*timePoint/(1.-timePoint/15000.))
        return c1

    def costSt2(timePoint, S,cost):
        c1 = costNoSt(timePoint)+cost
        aTemp = (S/(1.+10.*timePoint/(15000.+timePoint)))**2
        c2 = aTemp/(1. + aTemp)
        all = np.asarray(c1+c2-c1*c2)
        if all <0 or all > 0.9 or c1 < 0: all = 1
        growth = (1. - all - 0.1)/((1- cost-0.1)/(1-cost))
        if growth < 0: growth = 0

        return growth

    import numpy as np
    all = np.zeros([100,101])
    i = -1
    for DS in np.linspace(0,10000,100):
        i += 1
        j = -1
        for AB in np.linspace(0,10,101):
            j += 1
            all[i,j] = costSt2(DS,AB,cost)
            #if all[i,j] == 0: all[i,j] = -1
    print all
    np.savetxt('allCostsSt_W'+str(cost)+'.txt',all)
    import pylab as plt

def mainW(cost):
    def costNoSt(timePoint):
        c1 = np.asarray(1./100000.*timePoint/(1.-timePoint/15000.))
        return c1

    def costSt2(timePoint, S,cost):
        c1 = costNoSt(timePoint)+cost
        aTemp = (S/(1.+10.*timePoint/(15000.+timePoint)))**2
        c2 = aTemp/(1. + aTemp)
        all = np.asarray(c1+c2-c1*c2)
        if all <0 or all > 0.9 or c1 < 0: all = 1
        growth = (1. - all - 0.1)/((1- cost-0.1)/(1-cost))
        if growth < 0: growth = -1

        return growth

    import numpy as np
    all = np.zeros([100,101])
    i = -1
    for DS in np.linspace(0,10000,100):
        i += 1
        j = -1
        for AB in np.linspace(0,10,101):
            j += 1
            all[i,j] = costSt2(DS,AB,cost)
            #if all[i,j] == 0: all[i,j] = -1
    print all
    np.savetxt('allCostsSt_S'+str(cost)+'.txt',all)
    import pylab as plt

for x in [0,0.01,0.03,0.1]:
    main(x)

"""
from studyResults import customaxis
for x in [0,0.01,0.03,0.1]:
    #main(x)
    import pylab as plt
    import numpy as np
    all = np.loadtxt('allCostsSt_W'+str(x)+'.txt')
    fig, ax = plt.subplots()
    fig.set_size_inches(10.,10.)
    im = ax.imshow(all, cmap='BrBG', interpolation='none',
                   vmin=0, vmax=1,origin='lower')
    customaxis(ax, c_left='k', c_bottom='k', c_right='none', c_top='none', lw=2, size=20, pad=8)
    print np.linspace(0,10000,11)

    plt.xticks(np.linspace(0,100,6),[int(_) for _ in np.linspace(0,10,6)])
    plt.yticks(np.linspace(0,100,6),[int(_) for _ in np.linspace(0,10000,6)])
    plt.xlabel('Stressor [mM]',fontsize=20)
    plt.ylabel('Protein [molecules]',fontsize=20)
    fig.colorbar(im)
    plt.savefig('fitnessLandscape'+str(x)+'.pdf', bbox_inches='tight' ,dpi=100)
    plt.show()
"""
