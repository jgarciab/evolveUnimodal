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
"""
from problemBench import *
from pylab import *
def costNoSt(timePoint):
    c1 = 1./100000.*timePoint/(1.-timePoint/15000.)
    return c1

def  costSt(timePoint, S):
    c1 = costNoSt(timePoint)
    c2 = (S/(1.+10.*timePoint/(15000.+timePoint)))**2/(1. + (S/(1.+10.*timePoint/(15000.+timePoint)))**2)
    return c1+c2-c1*c2
    
    
c = [1000,0]
timeReaction = 1E6
A = c[1]+c[0]*(1+np.sin(np.linspace(0,timeReaction,timeReaction)/(12*np.pi)))

a = np.zeros([timeReaction,1])
for i in range(1):
    [mean,temp] = progress(timeReaction,250,10000.,1.,1.,actSeries = A,flag10Genes=1)
    print np.mean(mean),mstats.mquantiles(mean,prob=[0.9999])[0]
    a[:,i] = mean
a2 = np.nanmean(a,1)
plot(A,color='black',linewidth=2)
plot(a2)
show()

costNoStressor = np.mean(costNoSt(a2))

best = mstats.mquantiles(a2,prob=[0.9999])[0]
maxS = 0
print best, np.mean(a2), np.std(a2)
while costSt(best,maxS)<0.95:
    maxS += 0.01
#print np.mean(a2), np.max(a2), costNoStressor, -maxS
print costNoStressor, -maxS
"""

"""
import numpy as np
from pylab import *
Evol = np.loadtxt('EvolveNoiseFixed_'+str(1)+'GenesFitn_sinFreq.txt')
Evol2 = np.loadtxt('EvolveNoise'+str(1)+'GenesFitn_sinFreq.txt')
hold(True)
plot(Evol[:,0],Evol[:,1],'.',markersize=12)
plot(Evol2[:,0],Evol2[:,1],'.',markersize=12)
show()
"""

"""
import numpy as np
f3 = open('it.dic')

lines = [line for line in f3 if line.strip()]
lines.sort()
file = open("it2.dic", "w")
for i in lines:
    file.write(i)
file.close()

f3.close()
f3 = open('it2.dic')
a3 = f3.readline().split()
f2 = open('2.txt')
a2 = f2.readline().split()
f1 = open('1.csv')
count = 0
lines = [line for line in f1 if line.strip()]
f1.close()
lines.sort()

arrayStuff = np.zeros([5000,5000])
countW = np.zeros(5000)
print countW[1550]
arrayStuff[:,:] = np.NaN
legends = []
j = 0
jold = -1
for line in lines:
    a =line.split('\t')
    try:
        a[0].decode('ascii')
        if (a[0] < 'zzzz' and a[0] > 'ZZZ'): countW[int(a[1])] += int(a[2])
        if a[0]< a2[0]:
            pass
            #print '1', a[0],a2[0]
        elif a[0] == a2[0]:
            while a[0] > a3[0]:
                a3 = f3.readline().split()
                if len(a3) == 0: a3 = ['zzzzz']

            if a[0] == a3[0]:
                pass
            else:
                print '2', j, a[0],a2[0]
                arrayStuff[j,a[1]] = a[2]
    #            print a[1], a[2]

                if jold != j: legends.append(a[0])
                jold = j
        else:
            #print '3', a[0],a2[0]
            while a2[0] < a[0]:
                a2 = f2.readline().split()
                if len(a2) == 0: a2 = ['zzzzz']
            j+=1


    except UnicodeDecodeError:
        pass
print legends


for i in range(5000):
    arrayStuff[i,~np.isnan(arrayStuff[i,:])] = arrayStuff[i,~np.isnan(arrayStuff[i,:])]/countW[~np.isnan(arrayStuff[i,:])]


np.savetxt('3c.txt',arrayStuff)
file = open("legendsc.csv", "w")
for i in legends:
    file.write(i)
    file.write("\n")
file.close()


legends= []
file = open("legendsc.csv", "r")
for line in file:
    legends.append(line)
file.close()
print legends



arrayStuff = np.loadtxt('3c.txt')
from pylab import *
count = 0
pito = np.zeros(5000)

finalarray = np.zeros([512,1000])
t = 0
finalarray[:-2,t] = np.arange(1500,2010)
finalarray[-2:,t] = np.NaN
for i in range(5000):

    if ~all(np.isnan(arrayStuff[i,:])):
        print i
        t += 1
        finalarray[:-2,t] = arrayStuff[i,1500:2010]
        finalarray[-2,t] = np.nanmean(arrayStuff[i,:])
        finalarray[-1,t] = np.nanmax(arrayStuff[i,:])

        figure(1)
        plot(np.arange(5000)[arrayStuff[i,:] != np.NaN],arrayStuff[i,:])
        figure(2)
        plot(np.arange(5000)[~np.isnan(arrayStuff[i,:])],(arrayStuff[i,~np.isnan(arrayStuff[i,:])]-np.mean(arrayStuff[i,~np.isnan(arrayStuff[i,:])]))/np.std(arrayStuff[i,~np.isnan(arrayStuff[i,:])]))

        arrayStuff[i,~np.isnan(arrayStuff[i,:])] = (arrayStuff[i,~np.isnan(arrayStuff[i,:])]-np.mean(arrayStuff[i,~np.isnan(arrayStuff[i,:])]))/np.std(arrayStuff[i,~np.isnan(arrayStuff[i,:])])



        count+=1

savetxt('melc.csv',np.transpose(finalarray[:,:t+1]),delimiter=',')



figure(1)
xlabel('Time')
ylabel('Count')
#legend(legends,loc='center left', bbox_to_anchor=(1, 0.5))
savefig('11c.pdf', bbox_inches='tight' ,dpi=100)
figure(2)
xlabel('Time')
ylabel('Normalizzed count')
#legend(legends,loc='center left', bbox_to_anchor=(1, 0.5))
savefig('22c.pdf', bbox_inches='tight' ,dpi=100)

figure(3)
xlabel('Time')
ylabel('Count')
plot(np.arange(5000), nanmean(arrayStuff,0))
#legend(legends,loc='center left', bbox_to_anchor=(1, 0.5))
savefig('33c.pdf', bbox_inches='tight' ,dpi=100)

from scipy.stats import mstats
t = nonzero(pito > mstats.mquantiles(pito[pito!=0],prob=[0.95]))[0]
print t
print len(t), len(pito)
for cosa in t:
    print legends[cosa],
"""
def main(cost):
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

for x in y:
    main(x)

"""
from studyResults import customaxis
for x in [0]:#,0.01,0.03,0.1]:
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
