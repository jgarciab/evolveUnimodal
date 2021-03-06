import numpy as np
import pandas as pd
import pylab as plt

from studyResults import customaxis, ratio

# These are the "Tableau 20" colors as RGB.
cols = [(3,43,122),(31,119,180),(174,199,232),(255,186,120),(177,3,24),(31, 119, 180),
        (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229),
             (1, 119, 180), (1, 199, 232), (1, 127, 14), (1, 187, 120),
             (1, 160, 44), (1, 223, 138), (1, 39, 40), (1, 152, 150),
             (1, 103, 189), (1, 176, 213), (1, 86, 75), (1, 156, 148),
             (1, 119, 194), (1, 182, 210), (1, 127, 127), (1, 199, 199),
             (1, 189, 34), (1, 219, 141), (1, 190, 207), (1, 218, 229)]


# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.
for i in range(len(cols)):
    r, g, b = cols[i]
    cols[i] = (r / 255., g / 255., b / 255.)


def extractInfoFromFile(fnameF,noise,dim,candFlag = False):
    from scipy.stats import lognorm,gamma
    df = np.asarray(pd.read_csv(fnameF))

    maxI = len(df[:,1])-1

    cand = np.mean(df[-5:-1,3:-1],0)
    fit = df[-1,1]
    sense = df[-1,-1]

    if noise > 0:
        for line in range(10, len(df[:,1])):
            fit = 0.8*fit + 0.2*df[line,1]
            sense = 0.8*sense + 0.2*df[line,-1]
            cand = 0.9*cand + 0.1*df[line,3:-1]
            #print df[line,3:-1]

    if candFlag: return cand

    fit *= -1
    th = 100

    if dim == 100: d = cand/np.sum(cand)
    else:
        x = np.linspace(0.01, 10000., num=100) # values for x-axis
        d = np.zeros(100)
        w = 0
        for jj in range(0,len(cand)-1,3):
            d += cand[jj]*gamma.cdf(x, cand[jj+1], loc=0, scale=cand[jj+2]) # probability distribution
            w += cand[jj]
        d = np.diff(np.concatenate([[0],d]))
        d = d/w


    return [d,cand],sense,fit

def main2_4(tit,dim=1,noise=0,numE2=2):
    """
    version >4
    """
    

    x1 = [1,3,10,30,100]

    if dim == 100: end = "2STDobs.txt"
    else: end = "2obs.txt"
    numS = 5
    numE = 5
    allFit = np.zeros([numS,numE])
    allSk = np.zeros([numS,numE])
    allKut = np.zeros([numS,numE])
    allMe = np.zeros([numS,numE])
    allFa = np.zeros([numS,numE])
    allSt = np.zeros([numS,numE])
    allSense = np.zeros([numS,numE])
    propLowHigh = np.zeros([numS,numE])
    MeanAll = np.zeros([numS,numE,100])
    i = -1

    for a in range(numS) :
        i += 1
        j = -1
        for b in np.arange(numE2):
            j += 1
            #open("./dataDE/"+str(noise)+args+str(dim)+str(suddenness)+str(numChanges)+"obs.txt","w")
            fnameF = "./dataDE/"+str(noise)+tit+str(dim)+str(b)+str(x1[a])+end
            [cand,gammaP],sense,fit = extractInfoFromFile(fnameF,noise,dim)

            Mean = np.nansum(np.linspace(0,10000,100)*cand)
            Variance = np.nansum(cand*(np.linspace(0,10000,100) - Mean)**2)
            allSense[i,j] = np.round(sense)
            allFit[i,j] = fit
            MeanAll[i,j,:] = cand



            allMe[i,j] += Mean
            allFa[i,j] += np.sqrt(Variance/Mean)
            allSt[i,j] += np.sqrt(Variance)

            propLowHigh[i,j] = np.sum(MeanAll[i,j,40:])*100.
            if 0:#propLowHigh[i,j] == 0:
                plt.plot(cand)
                plt.title(str(a)+str(b))
                print tit
                plt.show()


    np.savetxt("./dataDE/EN_Oct19_allFit_v4"+str(noise)+str(dim)+tit+".txt",allFit)
    allFit[np.isnan(allFit)] = 0

    with file('./dataDE/EN_Oct19_MeanAll_v4'+str(noise)+str(dim)+tit+'.txt', 'w') as outfile:
        outfile.write('# Array shape: {0}\n'.format(MeanAll.shape))
        for data_slice in MeanAll:
            np.savetxt(outfile, data_slice, fmt='%-7.2e')
            outfile.write('# New slice\n')

    allMe[np.isnan(allMe)] = 0
    allFa[np.isnan(allFa)] = 0
    allSense[np.isnan(allSense)] = 0

    propLowHigh[np.isnan(propLowHigh)] = 0


    np.savetxt("./dataDE/EN_Oct19_allMe_v4"+str(noise)+str(dim)+tit+".txt",allMe)
    np.savetxt("./dataDE/EN_Oct19_allFa_v4"+str(noise)+str(dim)+tit+".txt",allFa)
    np.savetxt("./dataDE/EN_Oct19_allSense_v4"+str(noise)+str(dim)+tit+".txt",allSense)
    np.savetxt("./dataDE/EN_Oct19_propLowHigh_v4"+str(noise)+str(dim)+tit+".txt",propLowHigh)

def printAB(tit,dim=1,noise=0):
    try:
        x1 = [1,3,10,30,100]
        i = -1
        if dim == 100: end = "STDobs.txt"
        else: end = "obs.txt"
        a1,a2,a3,a4 = 0,0,0,0
        for a in range(5) :
            i += 1
            j = -1
            for b in np.arange(2):
                j += 1
                #open("./dataDE/"+str(noise)+args+str(dim)+str(suddenness)+str(numChanges)+"obs.txt","w")
                fnameF = "./dataDE/"+str(noise)+tit+str(dim)+str(b)+str(x1[a])+end

                df = np.asarray(pd.read_csv(fnameF))
                maxI = len(df[:,1])-1

                cand = df[-1,3:-1]
                fit = df[-1,1]
                sense = df[-1,-1]

                for line in range(10, len(df[:,1])):
                    fit = 0.8*fit + 0.2*df[line,1]
                    sense = 0.8*sense + 0.2*df[line,-1]
                    cand = 0.8*cand + 0.2*df[line,3:-1]

                a1 += fit
                a2 += sense
                a3 += cand[1]
                a4 += cand[2]
        print "%f, %f, %f, %f" %(fit,sense,cand[1],cand[2])
    except: pass

def plotImshow(numS,numE,MeanAll,allFit,allMe,allFa,allSense,propLowHigh,tit,limits=True):
    i = -1
    fig2 = plt.figure(2)
    fig2.set_size_inches(2.*numS,2.*numE)

    #10 Env for a in np.logspace(0.0,0.4,5)-1:#np.logspace(0.000,0.5,6)-1:

    for a in range(numS) :
        i += 1
        j = -1
        #10 Env for b in np.arange(5):#np.arange(0,7,1):#np.arange(1,10,1):
        for b in np.arange(numE):
            j += 1

            x =  numE*((numS-1) -i)+j+1
            #print x, i, j
            ax = fig2.add_subplot(numS,numE, x)

            ax.plot(np.linspace(50,10000,100),MeanAll[i,j,:],linewidth=3,color='orange')

            if i == 0:
                if j == 0:
                    ax.set_ylabel('Frequency',fontsize=16)
                    ax.set_xlabel('Protein \nLevel',fontsize=16)
                else:
                    ax.set_xlabel('Protein \nLevel',fontsize=16)
                    ax.axes.get_yaxis().set_visible(False)
            elif j == 0:
                ax.set_ylabel('Frequency',fontsize=16)
                ax.axes.get_xaxis().set_visible(False)
            else:
                ax.axes.get_xaxis().set_visible(False)
                ax.axes.get_yaxis().set_visible(False)
            frame = plt.gca()
            frame.axes.get_xaxis().set_ticks([])
            frame.axes.get_yaxis().set_ticks([])
            #ax.set_ylim((0,0.5))
            ax.set_xlim((0,10050))

            ax.set_xlim((-200,10000))
            plt.subplots_adjust(hspace = .001, wspace=0.001)



    tit = tit.replace('_','-')
    plt.suptitle(tit,fontsize=24)
    plt.savefig("./imagesDE/"+tit+'Hist.pdf', bbox_inches='tight' ,dpi=100)
    #plt.show()
    plt.clf()

def main2_4_Plot(tit,dim=1,noise=0,numE2=2):
    """
    version >4
    """
    x1 = [1,3,10,30,100]


    numS = 5
    numE = 5
    #print np.shape(np.loadtxt("./data/EN_Oct19_meanAll_v4"+tit+".txt"))
    MeanAll = np.loadtxt("./dataDE/EN_Oct19_MeanAll_v4"+str(noise)+str(dim)+tit+".txt").reshape((numS,numE,100))
    allFit = np.loadtxt("./dataDE/EN_Oct19_allFit_v4"+str(noise)+str(dim)+tit+".txt")
    allMe = np.loadtxt("./dataDE/EN_Oct19_allMe_v4"+str(noise)+str(dim)+tit+".txt")
    allFa = np.loadtxt("./dataDE/EN_Oct19_allFa_v4"+str(noise)+str(dim)+tit+".txt")
    allSense = np.loadtxt("./dataDE/EN_Oct19_allSense_v4"+str(noise)+str(dim)+tit+".txt")
    propLowHigh = np.loadtxt("./dataDE/EN_Oct19_propLowHigh_v4"+str(noise)+str(dim)+tit+".txt")


    plotImshow(numS,numE2,MeanAll,allFit,allMe,allFa,allSense,propLowHigh,str(noise)+str(dim)+tit)


def main2_4_PlotTogether(tit,dim=1,noise=0,numE2=2):
    """
    version >4
    """
    x1 = [1,3,10,30,100]

    numS = 5
    numE = 2



    fig2 = plt.figure(2)
    fig2.set_size_inches(2.,2.*5)
    jjj = 0
    #10 Env for a in np.logspace(0.0,0.4,5)-1:#np.logspace(0.000,0.5,6)-1:
    for d in dim:
        i = -1
        jjj += 1
        MeanAll = np.loadtxt("./dataDE/EN_Oct19_MeanAll_v4"+str(noise)+str(d)+tit+".txt").reshape((numS,5,100))
        for a in range(numS) :
            i += 1
            j = -1
            #10 Env for b in np.arange(5):#np.arange(0,7,1):#np.arange(1,10,1):
            for b in np.arange(numE):
                j += 1

                x =  numE*((numS-1) -i)+j+1
                #print x, i, j
                ax = fig2.add_subplot(numS,numE, x)

                ax.plot(np.linspace(50,10000,100),MeanAll[i,j,:],linewidth=3,color=cols[jjj])

                if i == 0:
                    if j == 0:
                        ax.set_ylabel('Frequency',fontsize=16)
                        ax.set_xlabel('Protein \nLevel',fontsize=16)
                    else:
                        ax.set_xlabel('Protein \nLevel',fontsize=16)
                        ax.axes.get_yaxis().set_visible(False)
                elif j == 0:
                    ax.set_ylabel('Frequency',fontsize=16)
                    ax.axes.get_xaxis().set_visible(False)
                else:
                    ax.axes.get_xaxis().set_visible(False)
                    ax.axes.get_yaxis().set_visible(False)
                frame = plt.gca()
                frame.axes.get_xaxis().set_ticks([])
                frame.axes.get_yaxis().set_ticks([])
                ax.set_ylim((0,0.05))
                ax.set_xlim((0,10050))

                ax.set_xlim((-200,10000))
                #plt.yscale('log')
                plt.subplots_adjust(hspace = .001, wspace=0.001)
                #plt.xscale('log')


    tit = tit.replace('_','-')
    plt.suptitle(tit,fontsize=24)

    plt.savefig("./imagesDE/"+tit+str(noise)+'Hist.pdf', bbox_inches='tight' ,dpi=100)

    #plt.show()
    plt.clf()

def figStrategiesRegion(titles,cost=1,dim=1,noise=0):
    numE = 5
    matrixToPlot = np.zeros((len(titles),numE))

    i = -1
    for tit in titles:
        i += 1
        propLowHigh = np.loadtxt("./dataDE/EN_Oct19_propLowHigh_v4"+str(noise)+str(dim)+tit+".txt")

        #matrixToPlot[i,:] = ((propLowHigh[:,cost]))#
        matrixToPlot[i,:] =  np.log10((propLowHigh[:,cost])/(100.-propLowHigh[:,cost]))

        """
        x = []
        for cosa in ['1','3','10','30','100']:
            fnameF = "./dataDE/"+str(noise)+tit+str(dim)+str(-1*(cost-1))+cosa+"0obs.txt"
            d = extractInfoFromFile(fnameF,noise,dim,candFlag = True)
            print d
            x.append(d[0]/d[3])

        print tit, x
        matrixToPlot[i,:] = np.asarray(x)
        """
    print matrixToPlot
    fig = plt.figure(1)
    fig.set_size_inches(7./5*2.,2.)

    ax = plt.subplot(1,1,1)
    #plt.imshow(matrixToPlot,interpolation='none',cmap='Blues',vmin=-20,vmax=100)
    if 1:#cost == 1:
        plt.imshow(matrixToPlot,interpolation='none',cmap='Purples',vmin=-5,vmax=3)#,vmin=-3,vmax=2)#,vmin=-10,vmax=100)
    else:
        plt.imshow(matrixToPlot,interpolation='none',cmap='Purples',vmin=-3,vmax=3)#,vmin=-3,vmax=2)#,vmin=-10,vmax=100)
    plt.xlabel('Environmental Transition Rate',fontsize=10)
    plt.ylabel('Symmetry',fontsize=10)
    plt.xticks(np.arange(5),['1','3','10','30','100'],rotation=90)



    plt.yticks(np.arange(len(titles)),["%.0e" %(float(ratio(title)[ratio(title).find(":")+1:])/float(ratio(title)[:ratio(title).find(":")])) for title in titles])
    plt.title('Proportion of cells in phenotype high',fontsize=10)
    plt.colorbar(orientation='vertical')
    customaxis(ax, c_left='none', c_bottom='none', c_right='none', c_top='none', lw=2, size=10, pad=8)
    plt.savefig("./imagesDE/imshow"+title[0]+str(cost)+str(noise)+str(dim)+tit[-2:]+'.pdf', bbox_inches='tight' ,dpi=100)
    plt.clf()
    #plt.show()

def figStrategiesPropTo(titles,cost=1,scale='log',dim=1,noise=0):
    numE = 5
    path = './data/'

    fig = plt.figure(1)
    fig.set_size_inches(3.3,2.5)
    ax = plt.subplot(1,1,1)
    minX = 100

    i = -1
    for tit in titles:
        i += 1
        propLowHigh = np.loadtxt("./dataDE/EN_Oct19_propLowHigh_v4"+str(noise)+str(dim)+tit+".txt")

        ax.plot([1,3,10,30,100],propLowHigh[:,cost],marker='.',ms=10,label=ratio(tit),linewidth=1.5,color = cols[i])
        t = propLowHigh[:,cost]
        minX = min(min(t[t>0]),minX)

    plt.xscale('log')
    plt.xlabel('Environmental Transition Rate',fontsize=10)
    plt.ylabel('Percentage of cells with \n high protein expression',fontsize=10)
    plt.yscale(scale)
    plt.ylim(minX/10.,120)
    plt.xlim(0.9,120)
    plt.legend(loc=4,frameon=0,fontsize=10,ncol=2)
    customaxis(ax, c_left='k', c_bottom='k', c_right='none', c_top='none', lw=1, size=10, pad=8)

    #plt.show()


    plt.show()
    plt.savefig("./imagesDE/lines"+tit[0]+str(cost)+str(noise)+str(dim)+tit[-2:]+'.pdf', bbox_inches='tight' ,dpi=100)
    plt.clf()

def figStrategiesAsymmetryo(titles,cost=1,scale='log',dim=1,noise=0):
    numE = 5
    path = './dataDE/'

    fig = plt.figure(1)
    fig.set_size_inches(0.8*3.3,0.8*2.5)
    ax = fig.add_subplot(1,1,1)



    x0 = []
    x1 = []
    x2 = []
    x3 = []
    x4 = []
    y0 = []
    y1 = []
    y2 = []
    y3 = []
    y4 = []
    i = -1
    xaxis = []
    xaxis2 = []


    for tit in titles:
        i += 1
        propLowHigh = (np.loadtxt("./dataDE/EN_Oct19_propLowHigh_v4"+str(noise)+str(dim)+tit+".txt"))
        print propLowHigh

        x0.append(propLowHigh[0,0])
        x1.append(propLowHigh[1,0])
        x2.append(propLowHigh[2,0])
        x3.append(propLowHigh[3,0])
        x4.append(propLowHigh[4,0])
        y0.append(propLowHigh[0,1])
        y1.append(propLowHigh[1,1])
        y2.append(propLowHigh[2,1])
        y3.append(propLowHigh[3,1])
        y4.append(propLowHigh[4,1])

        env = ratio(tit)
        if tit[0] == "0":
            xaxis.append(env)
        else:

            first = env[:env.find(":")]
            second = env[env.find(":")+1:]
            print float(env[:env.find(":")])/float(env[env.find(":")+1:])
            xaxis.append(float(env[env.find(":")+1:])/float(env[:env.find(":")]))
            xaxis2.append(float(env[:env.find(":")]))

    xaxis = np.asarray(xaxis)
    xaxis2 = np.asarray(xaxis2)
    """
    ax.plot(xaxis,x0,color=(214./255, 39./255, 40./255),marker='.',label='No Sensing, Env 1')
    ax.plot(xaxis,x1,color=(214./255, 39./255, 40./255),marker='.',label='No Sensing, Env 3')
    ax.plot(xaxis,x2,color=(214./255, 39./255, 40./255),marker='.',label='No Sensing, Env 10')
    ax.plot(xaxis,x3,color=(214./255, 39./255, 40./255),marker='.',label='No Sensing, Env 30')
    ax.plot(xaxis,x4,color=(214./255, 39./255, 40./255),marker='.',label='No Sensing, Env 100')

    ax.plot((xaxis),y0,color=(31./255, 119./255, 180./255),alpha=1,marker='.',label='Sensing, Env 1')
    ax.plot((xaxis),y1,color=(31./255, 119./255, 180./255),alpha=1,marker='.',label='Sensing, Env 3')
    ax.plot((xaxis),y2,color=(31./255, 119./255, 180./255),alpha=1,marker='.',label='Sensing, Env 10')
    ax.plot((xaxis),y3,color=(31./255, 119./255, 180./255),alpha=1,marker='.',label='Sensing, Env 30')
    ax.plot((xaxis),y4,color=(31./255, 119./255, 180./255),alpha=1.,marker='.',label='Sensing, Env 100')
    """
    ax.plot(1./(xaxis2/1),x0,color=(214./255, 39./255, 40./255),marker='.',label='No Sensing, Env 1')
    ax.plot(1./(xaxis2/3),x1,color=(214./255, 39./255, 40./255),marker='.',label='No Sensing, Env 3')
    ax.plot(1./(xaxis2/10),x2,color=(214./255, 39./255, 40./255),marker='.',label='No Sensing, Env 10')
    ax.plot(1./(xaxis2/30),x3,color=(214./255, 39./255, 40./255),marker='.',label='No Sensing, Env 30')
    ax.plot(1./(xaxis2/100),x4,color=(214./255, 39./255, 40./255),marker='.',label='No Sensing, Env 100')

    ax.plot((1./(xaxis2/1)),y0,color=(31./255, 119./255, 180./255),alpha=1,marker='.',label='Sensing, Env 1')
    ax.plot((1./(xaxis2/3)),y1,color=(31./255, 119./255, 180./255),alpha=1,marker='.',label='Sensing, Env 3')
    ax.plot((1./(xaxis2/10)),y2,color=(31./255, 119./255, 180./255),alpha=1,marker='.',label='Sensing, Env 10')
    ax.plot((1./(xaxis2/30)),y3,color=(31./255, 119./255, 180./255),alpha=1,marker='.',label='Sensing, Env 30')
    ax.plot((1./(xaxis2/100)),y4,color=(31./255, 119./255, 180./255),alpha=1.,marker='.',label='Sensing, Env 100')

    print y0
    """
    for i,j,z in zip(xaxis,y0,np.ones(len(y0))):
        print "%d  %d %d" %(i*1000,j*1000,z)
    for i,j,z in zip(xaxis,y1,2+np.ones(len(y0))):
        print "%d  %d %d" %(i*1000,j*1000,z)
    for i,j,z in zip(xaxis,y2,9+np.ones(len(y0)))[:-1]:
        print "%d  %d %d" %(i*1000,j*1000,z)
    for i,j,z in zip(xaxis,y3,29+np.ones(len(y0)))[:-1]:
        print "%d  %d %d" %(i*1000,j*1000,z)
    """
    ax.set_xlabel('Ratio High/Low Stress',fontsize=10)
    ax.set_xlabel('Time High Stress',fontsize=10)
    ax.set_ylabel('Percentage of cells with \n high protein expression',fontsize=10)

    ax.set_yscale(scale)
    ax.set_xscale(scale)
    #plt.ylim(1E-3,120)
    #plt.xlim(0.    ax2.plot(1./(xaxis*1),y0,color='blue',alpha=1,marker='.',label='Sensing, Env 1')
    #ax.legend(loc=4,frameon=0,fontsize=10,ncol=2)
    customaxis(ax, c_left='k', c_bottom='k', c_right='none', c_top='none', lw=1, size=10, pad=8)

    #plt.show()
    plt.savefig("./imagesDE/linesAsymmetry"+tit[0]+str(cost)+str(noise)+str(dim)+tit[-2:]+'.pdf', bbox_inches='tight' ,dpi=100)
    plt.clf()

def figStrategiesEnvironmentEffect(titles,noises = [0],cost= 0,dims = [2,1]):
    numE = 5
    path = './dataDE/'

    fig = plt.figure(1)
    fig.set_size_inches(3.3,2.5)
    ax = fig.add_subplot(1,1,1)

    if cost == 0:
        titles = ["2Env_NN_PEEHSS"]
        par = 2
    else:
        titles = ["2Env_NN_PEVHSS"]
        par = 4


    i = -1
    for noise in noises:

        for tit in titles:
            i += 1
            allFit = 100.*(np.loadtxt("./dataDE/EN_Oct19_allFit_v4"+str(noise)+str(dims[0])+tit+".txt")-np.loadtxt("./dataDE/EN_Oct19_allFit_v4"+str(noise)+str(dims[1])+tit+".txt"))#/np.loadtxt("./data/EN_Oct19_allFit_v4"+tit+".txt")
            propLowHigh = np.loadtxt("./dataDE/EN_Oct19_propLowHigh_v4"+str(noise)+str(dims[0])+tit+".txt")[par,cost]
            allFit = allFit[par,cost]
            ax.plot(propLowHigh.reshape(-1),allFit.reshape(-1),marker='.',ms=15,label=ratio(tit)+"N"+str(noise),color=cols[i])


    plt.ylim(-.1,2)
    plt.xlim(-5,105)
    ax.set_xlabel('Proportion of cells in High Stress')
    ax.set_ylabel(r'$\Delta Fitness$')
    plt.legend(loc=1,frameon=0,fontsize=10,ncol=2)
    customaxis(ax, c_left='k', c_bottom='k', c_right='none', c_top='none', lw=1, size=10, pad=8)
    #savefig('tempz.pdf', bbox_inches='tight' ,dpi=100)
    plt.savefig("./imagesDE/deltaFitness"+tit[0]+str(cost)+str(noise)+tit[-2:]+'.pdf', bbox_inches='tight' ,dpi=100)
    plt.clf()

def figStrategiesDFitnessSensing(titles,cost= [1,0],dims = 100):
    numE = 5
    path = './dataDE/'

    fig = plt.figure(1)
    fig.set_size_inches(3.3,2.5)
    ax = fig.add_subplot(1,1,1)
    i = -1

    for tit in titles:
        i += 1
        allFit = 100.*(np.loadtxt("./dataDE/EN_Oct19_allFit_v4"+str(0)+str(dims)+tit+".txt"))
        propLowHigh = np.loadtxt("./dataDE/EN_Oct19_propLowHigh_v4"+str(0)+str(dims)+tit+".txt")[:,1]

        allFit = allFit[:,1]-allFit[:,0]

        ax.plot(propLowHigh.reshape(-1),allFit.reshape(-1),marker='.',ms=15,label=ratio(tit)+"N"+str(0),color=cols[i])


    #plt.ylim(-.1,2)
    plt.xlim(-5,105)
    ax.set_xlabel('Proportion of cells with \n high protein expression')
    ax.set_ylabel(r'$\Delta Fitness$')
    plt.legend(loc=1,frameon=0,fontsize=10,ncol=2)
    customaxis(ax, c_left='k', c_bottom='k', c_right='none', c_top='none', lw=1, size=10, pad=8)
    #savefig('tempz.pdf', bbox_inches='tight' ,dpi=100)
    plt.savefig("./imagesDE/deltaFitnessSensing"+tit[0]+str(cost)+str(0)+tit[-2:]+'.pdf', bbox_inches='tight' ,dpi=100)
    plt.clf()

    numE = 5
    path = './dataDE/'

    fig = plt.figure(1)
    fig.set_size_inches(0.8*4,0.8*2.5)
    ax = fig.add_subplot(1,1,1)


    y0 = []
    y1 = []
    y2 = []
    y3 = []
    y4 = []
    i = -1
    xaxis = []
    xaxis2 = []


    for tit in titles:
        i += 1
        propLowHigh = 100.*(np.loadtxt("./dataDE/EN_Oct19_allFit_v4"+str(0)+str(dims)+tit+".txt"))

        propLowHigh[:,1] = propLowHigh[:,1]-propLowHigh[:,0]

        y0.append(propLowHigh[0,1])
        y1.append(propLowHigh[1,1])
        y2.append(propLowHigh[2,1])
        y3.append(propLowHigh[3,1])
        y4.append(propLowHigh[4,1])

        env = ratio(tit)
        if tit[0] == "0":
            xaxis.append(env)
        else:

            first = env[:env.find(":")]
            second = env[env.find(":")+1:]
            print float(env[:env.find(":")])/float(env[env.find(":")+1:])
            xaxis.append(float(env[env.find(":")+1:])/float(env[:env.find(":")]))
            xaxis2.append(float(env[:env.find(":")]))

    xaxis = np.asarray(xaxis)
    xaxis2 = np.asarray(xaxis2)

    ax.plot(xaxis,y0,color=(31./255, 119./255, 180./255),alpha=1,marker='.',label='Sensing, Env 1')
    ax.plot(xaxis,y1,color=(31./255, 119./255, 180./255),alpha=1,marker='.',label='Sensing, Env 3')
    ax.plot(xaxis,y2,color=(31./255, 119./255, 180./255),alpha=1,marker='.',label='Sensing, Env 10')
    ax.plot(xaxis,y3,color=(31./255, 119./255, 180./255),alpha=1,marker='.',label='Sensing, Env 30')
    ax.plot(xaxis,y4,color=(31./255, 119./255, 180./255),alpha=1.,marker='.',label='Sensing, Env 100')

    plt.ylim([0,12])
    ax.set_xlabel('Ratio High/Low Stress',fontsize=12)

    ax.set_ylabel(r'$\Delta Fitness$',fontsize=12)

    ax.set_yscale('linear')
    ax.set_xscale('log')

    customaxis(ax, c_left='k', c_bottom='k', c_right='none', c_top='none', lw=1, size=10, pad=8)

    #plt.show()
    plt.savefig("./imagesDE/linesAsymmetrxxxy"+tit[0]+str(cost)+str(0)+str(100)+tit[-2:]+'.pdf', bbox_inches='tight' ,dpi=100)
    plt.clf()

def figNoiseFitnessDifference(names,noises,cost=1,dims=[2,1]):
    numE = 5
    path = './dataDE/'

    fig = plt.figure(1)
    fig.set_size_inches(3.3,2.5)
    ax = fig.add_subplot(1,1,1)


    colsInd = -1
    for title in names:

        colsInd += 1
        i = -1
        structureFitness = np.zeros((len(noises),5))

        for noise in noises:
            i += 1
            print np.loadtxt("./dataDE/EN_Oct19_allFit_v4"+str(noise)+str(dims[0])+title+".txt")
            print np.loadtxt("./dataDE/EN_Oct19_allFit_v4"+str(noise)+str(dims[1])+title+".txt")
            allFit = 100.*(np.loadtxt("./dataDE/EN_Oct19_allFit_v4"+str(noise)+str(dims[0])+title+".txt")-
                           np.loadtxt("./dataDE/EN_Oct19_allFit_v4"+str(noise)+str(dims[1])+title+".txt"))#/np.loadtxt("./data/EN_Oct19_allFit_v4"+tit+".txt")
            structureFitness[i,:] = allFit[:,cost].reshape(-1)

        plt.title(title.replace("_","-"))
        plt.savefig("./imagesDE/tempTogether3Env"+title.replace("_","-")+str(cost)+title[-2:]+'.pdf', bbox_inches='tight' ,dpi=100)
        plt.clf()


        for j in range(5):
            if j ==  0:
                ax.plot(range(len(noises)),structureFitness[:,j],linewidth=2,color=cols[colsInd],label=ratio(title))
            else:
                ax.plot(range(len(noises)),structureFitness[:,j],linewidth=2,color=cols[colsInd])

    ax.yaxis.set_ticks_position('both')

    ax.vlines(range(len(noises)),0,2,linewidth=2,alpha=0.5)
    plt.ylim((-0.25,2.1))
    plt.xticks(range(len(noises)),noises)
    ax.set_xlabel('Noise level')
    ax.set_ylabel(r'$\Delta Fitness$')
    customaxis(ax, c_left='k', c_bottom='none', c_right='none', c_top='none', lw=2, size=10, pad=8)
    #plt.ylim(0,3)
    plt.legend(loc=1,frameon=0,fontsize=10,ncol=2)
    #savefig('tempw.pdf', bbox_inches='tight' ,dpi=100)
    plt.savefig("./imagesDE/deltaFitnessNoise"+title[0]+str(cost)+title[-2:]+'.pdf', bbox_inches='tight' ,dpi=100)
    plt.clf()

def figNoiseFitnessDifferenceOnlyOptimum(noises,cost=0,dims=[2,1]):
    numE = 5
    path = './dataDE/'
    import statsmodels.api as sm
    lowess = sm.nonparametric.lowess

    fig = plt.figure(1)
    fig.set_size_inches(3.5,2.5)
    ax = fig.add_subplot(1,1,1)

    noS = []
    S = []


    noS.append((10,"2Env_NN_PEAHSS"))
    noS.append((10,"2Env_NN_PEEHSS"))
    noS.append((10,"2Env_NN_PEVHSS"))
    noS.append((10,"2Env_NN_PEHHSS"))
    noS.append((10,"2Env_NN_PELHSS"))

    S.append((100,"2Env_NN_PEAHSS"))
    S.append((100,"2Env_NN_PEEHSS"))
    S.append((100,"2Env_NN_PEVHSS"))
    S.append((10,"2Env_NN_PEHHSS"))
    S.append((1,"2Env_NN_PELHSS"))


    allP = [noS,S]
    param = allP[cost]

    colsInd = -1
    for freq,title in param:
        colsInd += 1
        i = -1
        structureFitness = np.zeros((len(noises),6))
        structureFitness1 = np.zeros((len(noises),1))

        for noise in noises:
            i += 1
            fit0 = []
            plt.subplot(10,1,i+1)
            for iii in range(0,1):#[3,5]:
                fnameF = "./dataDE/"+str(noise)+title+str(dims[0])+str(cost)+str(freq)+str(iii)+"obs.txt"
                [cand,gammaP1],sense,fit = extractInfoFromFile(fnameF,noise,dims[0])

                plt.plot(cand)
                fit0.append(fit)
            fit1 = []
            for iii in range(0,1):#:
                fnameF = "./dataDE/"+str(noise)+title+str(dims[1])+str(cost)+str(freq)+str(iii)+"obs.txt"
                [cand,gammaP2],sense,fit = extractInfoFromFile(fnameF,noise,dims[1])

                plt.plot(cand)
                fit1.append(fit)

            #print "2g: ", gammaP1, "1g: ", gammaP2

            fit0 = np.asarray(fit0)
            fit1 = np.asarray(fit1)

            print cost, title[-5:-2], noise, fit0-fit1

            structureFitness[i,:] = 100.*(fit0-fit1)
            structureFitness1[i] = 100.*(np.median(fit0)-np.median(fit1))

            #structureFitness[i] = 100.*np.max(fit0)
            #structureFitness1[i] = 100.*np.max(fit1)

        def smooth(y, box_pts):
            box = np.ones(box_pts)/box_pts
            y_smooth = np.convolve(y, box, mode='same')
            return y_smooth

        ##plt.plot(np.repeat(noises,1),structureFitness.reshape(-1))
        ##plt.plot(np.repeat(noises,1),structureFitness1.reshape(-1))
        #ax.plot(np.repeat(noises,6),structureFitness.reshape(-1),ms=10,marker='.',color=cols[colsInd],linewidth=1.5)
        ##ax.plot(np.repeat(noises,1),structureFitness1.reshape(-1),color=cols[colsInd]
        ##a = lowess(structureFitness1.reshape(-1),np.repeat(noises,1),frac = 0.6,delta=0.0)

        ##ax.plot(np.repeat(noises,6),structureFitness.reshape(-1),'.',color=cols[colsInd])
        ##ax.plot(a[:,0],a[:,1],linewidth=2,color=cols[colsInd])
        plt.title(title.replace("_","-"))
        plt.savefig("./imagesDE/tempTogetherNoise"+title.replace("_","-")+str(cost)+title[-2:]+'.pdf', bbox_inches='tight' ,dpi=100)
        plt.clf()

    ax.yaxis.set_ticks_position('both')

    #ax.vlines(range(len(noises)),0,2,linewidth=2,alpha=0.5)
    plt.ylim((-0.1,2.6))
    #plt.xticks(range(len(noises)),noises)
    ax.set_xlabel('Noise level')
    ax.set_ylabel('Benefit of multistability')
    customaxis(ax, c_left='k', c_bottom='none', c_right='none', c_top='none', lw=2, size=10, pad=8)
    #plt.ylim(0,3)
    #plt.legend(loc=1,frameon=0,fontsize=10,ncol=2)

    plt.savefig("./imagesDE/deltaF2itnessNoiseOnlOptimum"+title[0]+str(cost)+title[-2:]+'.pdf', bbox_inches='tight' ,dpi=100)

    plt.clf()

def figNoiseFitnessDifferenceOnlyOptimum3Env(noises,cost=0,dims=[2,1]):
    numE = 5
    path = './dataDE/'
    import statsmodels.api as sm
    lowess = sm.nonparametric.lowess

    fig = plt.figure(1)
    fig.set_size_inches(3.5,2.5)
    ax = fig.add_subplot(1,1,1)

    noS = []
    S = []


    noS.append((10,"3Env_0102_PEA"))
    noS.append((10,"3Env_0102_PEE"))
    noS.append((10,"3Env_0102_PEV"))
    noS.append((10,"3Env_0102_PEH"))
    noS.append((10,"3Env_0102_PEL"))

    S.append((100,"3Env_0102_PEA"))
    S.append((100,"3Env_0102_PEE"))
    S.append((100,"3Env_0102_PEV"))
    S.append((10,"3Env_0102_PEH"))
    S.append((1,"3Env_0102_PEL"))




    allP = [noS,S]
    param = allP[cost]

    colsInd = -1
    for freq,title in param:
        colsInd += 1
        i = -1
        structureFitness = np.zeros((len(["A","E","V","H","L"]),1))
        structureFitness1 = np.zeros((len(["A","E","V","H","L"]),1))

        for ending in ["A","E","V","H","L"][::-1]:
            i += 1
            plt.subplot(5,1,i+1)
            fit0 = []
            for iii in [0]:#"["",1,2,3,4,5]:
                fnameF = "./dataDE/"+str(0)+title+str(ending)+"HSS"+str(dims[0])+str(cost)+str(freq)+str(iii)+"obs.txt"
                [cand,gammaP],sense,fit = extractInfoFromFile(fnameF,0,dims[0])
                plt.plot(cand)
                fit0.append(fit)
            fit1 = []
            for iii in [0]:#["",1,2,3,4,5]:
                fnameF = "./dataDE/"+str(0)+title+str(ending)+"HSS"+str(dims[1])+str(cost)+str(freq)+str(iii)+"obs.txt"
                [cand,gammaP],sense,fit = extractInfoFromFile(fnameF,0,dims[1])
                plt.plot(cand)
                fit1.append(fit)

            fit0 = np.asarray(fit0)
            fit1 = np.asarray(fit1)
            print fit0,fit1
            structureFitness[i,:] = 100.*(fit0-fit1)
            #structureFitness1[i,:] = 100.*

        plt.title(title.replace("_","-"))
        plt.savefig("./imagesDE/tempTogether3Env"+title.replace("_","-")+str(cost)+title[-2:]+'.pdf', bbox_inches='tight' ,dpi=100)
        plt.clf()

        structure = [4,3,2,1,0][::-1]
        #plt.plot(np.repeat(len(noises),1),structureFitness.reshape(-1))
        #plt.plot(np.repeat(len(noises),1),structureFitness1.reshape(-1))
        #a = lowess(structureFitness.reshape(-1),np.repeat(structure,6),frac = 0.66,delta=0.0)
        ax.plot(np.repeat(structure,1),structureFitness.reshape(-1),ms=10,marker='.',color=cols[colsInd],linewidth=1.5)
        #ax.plot(a[:,0],a[:,1],linewidth=2,color=cols[colsInd])


    ax.yaxis.set_ticks_position('both')

    #ax.vlines(range(len(noises)),0,2,linewidth=2,alpha=0.5)
    plt.ylim((-0.1,2.6))
    plt.xticks(range(len(structure)),structure)
    ax.set_xlabel('Time in the Intermediate Environment')
    ax.set_ylabel('Benefit of Multistability')
    customaxis(ax, c_left='k', c_bottom='none', c_right='none', c_top='none', lw=2, size=10, pad=8)
    #plt.ylim(0,3)
    #plt.legend(loc=1,frameon=0,fontsize=10,ncol=2)

    plt.savefig("./imagesDE/deltaF3ENVitnessNoiseOnlOptimum"+title[0]+str(cost)+title[-2:]+'.pdf', bbox_inches='tight' ,dpi=100)

    plt.clf()

def fig3Strategies(names,cost=1,dims=[2,1]):
    numE = 5
    path = './dataDE/'

    fig = plt.figure(1)
    fig.set_size_inches(3.3,2.5)
    ax = fig.add_subplot(1,1,1)

    colsInd = -1
    for n in ["A","E","V","H","L"][::-1]:

        n2 = [_ for _ in names if "3Env_Noise_PE"+n in _][::-1]
        colsInd += 1
        i = -1
        structureFitness = np.zeros((len(["A","E","V","H","L"]),5))

        for title in n2:
            i += 1
            allFit = 100.*(np.loadtxt("./dataDE/EN_Oct19_allFit_v4"+str(0)+str(dims[0])+title+".txt")-
                           np.loadtxt("./dataDE/EN_Oct19_allFit_v4"+str(0)+str(dims[1])+title+".txt"))#/np.loadtxt("./data/EN_Oct19_allFit_v4"+tit+".txt")


            structureFitness[i,:] = allFit[:,cost].reshape(-1)

        if n == "A": lab= 10000
        elif n == "E": lab= 1000
        elif n == "V": lab= 100
        elif n == "H": lab= 10
        elif n == "L": lab= 1

        for j in range(5):
            if j ==  0:
                ax.plot(range(5),structureFitness[:,j],linewidth=2,color=cols[len(["A","E","V","H","L"])-colsInd],label=str(lab))
            else:
                ax.plot(range(5),structureFitness[:,j],linewidth=2,color=cols[len(["A","E","V","H","L"])-colsInd])

    ax.yaxis.set_ticks_position('both')

    ax.vlines(range(5),0,2,linewidth=2,alpha=0.5)
    plt.xlim((0,3.1))
    plt.ylim((-0.25,2.1))

    plt.xticks(range(5),["1","10","100","1000","10000"])
    ax.set_xlabel('Second environment')
    ax.set_ylabel(r'$\Delta Fitness$')
    customaxis(ax, c_left='k', c_bottom='none', c_right='none', c_top='none', lw=2, size=10, pad=8)
    #plt.ylim(0,3)
    plt.legend(loc=1,frameon=0,fontsize=10,ncol=2)
    plt.savefig("./imagesDE/deltaFitnessSecondEnvir"+title[0]+str(cost)+str(0)+title[-2:]+'.pdf', bbox_inches='tight' ,dpi=100)
    plt.clf()

def figStrategiesRegionDeltaFitness(titles,cost=1,dims = [2,1]):
    numE = 5
    matrixToPlot = np.zeros((len(titles),numE))

    i = -1
    for title in titles:
        i += 1
        allFit = 100.*(np.loadtxt("./dataDE/EN_Oct19_allFit_v4"+str(0)+str(dims[0])+title+".txt")-
                       np.loadtxt("./dataDE/EN_Oct19_allFit_v4"+str(0)+str(dims[1])+title+".txt"))#/np.loadtxt("./data/EN_Oct19_allFit_v4"+tit+".txt")

        matrixToPlot[i,:] = (allFit[:,cost])

    print matrixToPlot
    fig = plt.figure(1)
    fig.set_size_inches(len(titles)/5.*2.,2.*2.)

    ax = fig.add_subplot(1,1,1)
    plt.imshow(matrixToPlot,interpolation='none',cmap='YlGn',vmin=-0.,vmax=1)
    #plt.imshow(matrixToPlot,interpolation='none',cmap='BuPu',vmin=-0.2,vmax=1.5)


    plt.xlabel('Environmental Transition Rate',fontsize=10)
    plt.ylabel('Symmetry',fontsize=10)
    plt.xticks(np.arange(5),['1','3','10','30','100'],rotation=90)


    print [float(ratio(title)[ratio(title).find(":")+1:])/float(ratio(title)[:ratio(title).find(":")]) for title in titles]
    plt.yticks(np.arange(len(titles)),[float(ratio(title)[ratio(title).find(":")+1:])/float(ratio(title)[:ratio(title).find(":")]) for title in titles])
    plt.title('Proportion of cells in phenotype high',fontsize=10)
    plt.colorbar(orientation='vertical')
    customaxis(ax, c_left='none', c_bottom='none', c_right='none', c_top='none', lw=2, size=10, pad=8)
    plt.savefig("./imagesDE/imshow3DeltaFit"+title[0]+str(cost)+str(0)+str(21)+title[-2:]+'.pdf', bbox_inches='tight' ,dpi=100)
    plt.clf()
    #plt.show()

def figGamma(names,dims=[2,1]):
    for cost in [0,1]:
        numE = 5
        path = './dataDE/'
        import statsmodels.api as sm
        lowess = sm.nonparametric.lowess

        fig = plt.figure(1)
        fig.set_size_inches(3.5,2.5)
        ax = fig.add_subplot(1,1,1)

        noS = []
        S = []

        fit0 = []
        fit1 = []
        colsInd = -1
        for title in names:
            freq = 100
            colsInd += 1
            i = -1


            for ending in [""][::-1]:
                i += 1
                for iii in [0]:#"["",1,2,3,4,5]:
                    fnameF = "./dataDE/"+str(0)+title+str(ending)+str(dims[0])+str(cost)+str(freq)+"obs.txt"
                    [cand,gammaP],sense,fit = extractInfoFromFile(fnameF,0,dims[0])
                    fit0.append(fit)

                for iii in [0]:#["",1,2,3,4,5]:
                    fnameF = "./dataDE/"+str(0)+title+str(ending)+str(dims[1])+str(cost)+str(freq)+"obs.txt"
                    [cand,gammaP],sense,fit = extractInfoFromFile(fnameF,0,dims[1])
                    fit1.append(fit)

        fit0 = np.asarray(fit0)
        fit1 = np.asarray(fit1)

        arrayDiffFitness = 100.*(fit0-fit1)


        st = np.array([0.001,0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
        st = st/1.5*11/10
        print len(st),len(arrayDiffFitness)
        #plt.plot(np.repeat(len(noises),1),structureFitness.reshape(-1))
        #plt.plot(np.repeat(len(noises),1),structureFitness1.reshape(-1))
        #a = lowess(structureFitness.reshape(-1),np.repeat(structure,6),frac = 0.66,delta=0.0)
        ax.plot(st,arrayDiffFitness.reshape(-1),ms=10,marker='.',color=cols[cost],linewidth=1.5,label=str(cost))
        #ax.plot(a[:,0],a[:,1],linewidth=2,color=cols[colsInd])


    ax.yaxis.set_ticks_position('both')

    #ax.vlines(range(len(noises)),0,2,linewidth=2,alpha=0.5)
    plt.ylim((-0.1,1))
    plt.xscale('log')
    plt.xlim((6E-4,1.))
    #plt.xticks(range(len(structure)),structure)
    ax.set_xlabel('Mean stress levels')
    ax.set_ylabel('Benefit of bistability')
    customaxis(ax, c_left='k', c_bottom='none', c_right='none', c_top='none', lw=2, size=10, pad=8)
    #plt.ylim(0,3)
    plt.legend(loc=1,frameon=0,fontsize=10,ncol=2)

    plt.savefig("./imagesDE/deltaFGamma"+title[0]+str(cost)+title[-2:]+'.pdf', bbox_inches='tight' ,dpi=100)
    plt.show()
    plt.clf()


def extractOptimum1Env():
    from scipy.stats import gamma
    gammaParameters = np.zeros((2,101))
    for stress in range(0,101):
        if stress < 10:
            s = "0"+str(stress)
        else:
            s = str(stress)
        title = "1Env_"+s+"SS"
        fnameF = "./dataDE/"+str(0)+title+"1010obs.txt"
        [d,cand],sense,fit = extractInfoFromFile(fnameF,0,1)

        x = np.linspace(0.01, 10000., num=100) # values for x-axis
        d = np.zeros(100)
        w = 0
        for jj in range(0,len(cand)-1,3):
            d += cand[jj]*gamma.pdf(x, cand[jj+1], loc=0, scale=cand[jj+2]) # probability distribution
            w += cand[jj]

        d /= np.sum(d)
        plt.plot(d)
        gammaParameters[:,stress] = cand[1:]
        print title, cand, fit
    plt.show()
    np.savetxt("gamma1EnvOptimum.txt",gammaParameters)


def main():

  """
  ## Fig. 2 
  names = ["2Env_NN_PEAHSS","2Env_NN_PEEHSS","2Env_NN_PEVHSS","2Env_NN_PEHHSS","2Env_NN_PELHSS"]
  possibleFactors = []
  for numChanges in  [1,3,10,30,100]:
      for sudden in range(2):
          for dim in [1,2]:
              for noise in [0]:
                  main2_4(name,dim=dim,noise=noise,numE2=1)
                  #main2_4_Plot(name,dim=dim,noise=noise,numE2=2)
              figStrategiesRegion(names,cost,dim,noise=0) #imshow
  figStrategiesDFitnessSensing(names,cost=[1,0],dims=2) #lognormal difference as propto in xcaxis
  """
  """
  ## Fig. 3 
  names = ["2Env_NN_PEAHSS","2Env_NN_PEEHSS","2Env_NN_PEVHSS","2Env_NN_PEHHSS","2Env_NN_PELHSS"]
  possibleFactors = []
  for numChanges in  [1,3,10,30,100]:
      for sudden in range(2):
          for dim in [1,2]:
              for noise in [0]:
                  main2_4(name,dim=dim,noise=noise,numE2=1)
                  #main2_4_Plot(name,dim=dim,noise=noise,numE2=2)
          figStrategiesRegionDeltaFitness(names,cost,dims=[2,1])
  """

  """
  ## Fig. 4 
  for cost in [0,1]:
      ## To plot distributions uncomment: subplot(...) and plot(cand) in the following files.	

      ## For 3 environments, difference in fitness as a function of the second environment. Only for optimum
      figNoiseFitnessDifferenceOnlyOptimum3Env([],cost=cost,dims=[2,1])

      ## Noise only for optimum
      figNoiseFitnessDifferenceOnlyOptimum([0,0.25, 0.5,0.75,1,1.5,2,3,4,5],cost=cost)
  """

  """
  ## Fig. S2 
  names = ["2Env_NN_PEAHSS","2Env_NN_PEEHSS","2Env_NN_PEVHSS","2Env_NN_PEHHSS","2Env_NN_PELHSS"]
  possibleFactors = []
  for numChanges in  [1,3,10,30,100]:
      for sudden in range(2):
          for dim in [100]:
              for noise in [0]:
                  main2_4(name,dim=dim,noise=noise,numE2=1)
                  #main2_4_Plot(name,dim=dim,noise=noise,numE2=2)
              figStrategiesRegion(names,cost,dim,noise=0) #imshow
  """

  """
  ## Fig. S3 
  names = ["2Env_NN_PEAHWS","2Env_NN_PEEHWS","2Env_NN_PEVHWS","2Env_NN_PEHHWS","2Env_NN_PELHWS"]
  possibleFactors = []
  for numChanges in  [1,3,10,30,100]:
      for sudden in range(2):
          for dim in [1,2]:
              for noise in [0]:
                  main2_4(name,dim=dim,noise=noise,numE2=1)
                  #main2_4_Plot(name,dim=dim,noise=noise,numE2=2)
              figStrategiesRegion(names,cost,dim,noise=0) #imshow
  figStrategiesDFitnessSensing(names,cost=[1,0],dims=2) #lognormal difference as propto in xcaxis
  """

  """
  ## Fig. S3 insets
  names = ["2Env_NN_PEIHWS","2Env_NN_PEtHWS","2Env_NN_PEjHWS","2Env_NN_PEkHWS","2Env_NN_PEsHWS","2Env_NN_PEmHWS", "2Env_NN_PEnHWS"]
  possibleFactors = []
  for numChanges in  [1,3,10,30,100]:
      for sudden in range(2):
          for dim in [1,2]:
              for noise in [0]:
                  main2_4(name,dim=dim,noise=noise,numE2=1)
                  #main2_4_Plot(name,dim=dim,noise=noise,numE2=2)
              figStrategiesRegion(names,cost,dim,noise=0) #imshow
  """
  """
   ## Fig. S5
   names = ["3Env_0102_PEHAHSS","3Env_0102_PEHEHSS","3Env_0102_PEHVHSS","3Env_0102_PEHHHSS","3Env_0102_PEHLHSS"]
   figGamma(names,dims=[2,1])
  """

  """
  ## OTHER STUFF.
    ## Extract data into something easily studied

    names =   ["3Env_Noise_PEAAHSS","3Env_Noise_PEAEHSS","3Env_Noise_PEAVHSS","3Env_Noise_PEAHHSS","3Env_Noise_PEALHSS",
                "3Env_Noise_PEEAHSS","3Env_Noise_PEEEHSS","3Env_Noise_PEEVHSS","3Env_Noise_PEEHHSS","3Env_Noise_PEELHSS",
                "3Env_Noise_PEVAHSS","3Env_Noise_PEVEHSS","3Env_Noise_PEVVHSS","3Env_Noise_PEVHHSS","3Env_Noise_PEVLHSS",
                "3Env_Noise_PEHAHSS","3Env_Noise_PEHEHSS","3Env_Noise_PEHVHSS","3Env_Noise_PEHHHSS","3Env_Noise_PEHLHSS",
                "3Env_Noise_PELAHSS","3Env_Noise_PELEHSS","3Env_Noise_PELVHSS","3Env_Noise_PELHHSS","3Env_Noise_PELLHSS"]

    names =   ["3Env_0102_PEAAHSS","3Env_0102_PEAEHSS","3Env_0102_PEAVHSS","3Env_0102_PEAHHSS","3Env_0102_PEALHSS",
                "3Env_0102_PEEAHSS","3Env_0102_PEEEHSS","3Env_0102_PEEVHSS","3Env_0102_PEEHHSS","3Env_0102_PEELHSS",
                "3Env_0102_PEVAHSS","3Env_0102_PEVEHSS","3Env_0102_PEVVHSS","3Env_0102_PEVHHSS","3Env_0102_PEVLHSS",
                "3Env_0102_PEHAHSS","3Env_0102_PEHEHSS","3Env_0102_PEHVHSS","3Env_0102_PEHHHSS","3Env_0102_PEHLHSS",
                "3Env_0102_PELAHSS","3Env_0102_PELEHSS","3Env_0102_PELVHSS","3Env_0102_PELHHSS","3Env_0102_PELLHSS"]

    #names = ["2Env_NN_PEAHSS","2Env_NN_PEEHSS","2Env_NN_PEVHSS","2Env_NN_PEHHSS","2Env_NN_PELHSS"]
    #names = sorted(["0Env_00010SS","0Env_0010SS","0Env_010SS","0Env_030SS","0Env_050SS","0Env_070SS","0Env_090SS","0Env_100SS","0Env_020SS","0Env_040SS","0Env_060SS","0Env_080SS"])
    #names = ["2Env_NN_PEAHWS","2Env_NN_PEEHWS","2Env_NN_PEVHWS","2Env_NN_PEHHWS","2Env_NN_PELHWS"]
    #names = ["2Env_NN_PEIHWS","2Env_NN_PEtHWS","2Env_NN_PEjHWS","2Env_NN_PEkHWS","2Env_NN_PEsHWS","2Env_NN_PEmHWS", "2Env_NN_PEnHWS"]
    #names = ["2Env_NN_PEAHSS","2Env_NN_PEEHSS","2Env_NN_PEVHSS","2Env_NN_PEHHSS","2Env_NN_PELHSS"]
    #names = ["2Env_NN_PEVHSS"]
    #names = ["2Env_NN_PEAHWS","2Env_NN_PEEHWS","2Env_NN_PEVHWS","2Env_NN_PEHHWS","2Env_NN_PELHWS"]


    """
    names = ["3Env_0102_PEHLHSS"]

    names = []
    for end in ["A","E","V","H","L"]:
        names.append("3Env_0102_PEE"+end+"HSS")
        names.append("3Env_0102_PEV"+end+"HSS")
    """

    #names = ["2Env_NN_PEAHSS","2Env_NN_PEEHSS","2Env_NN_PEVHSS","2Env_NN_PEHHSS","2Env_NN_PELHSS"]
    names = ["2Env_NN_PEIHWS","2Env_NN_PEtHWS","2Env_NN_PEjHWS","2Env_NN_PEkHWS","2Env_NN_PEsHWS","2Env_NN_PEmHWS", "2Env_NN_PEnHWS"]
    #names = ["2Env_NN_PEAHWS","2Env_NN_PEEHWS","2Env_NN_PEVHWS","2Env_NN_PEHHWS","2Env_NN_PELHWS"]
    #extractOptimum1Env()

    #names = ["2Env_NN_PEIHWS","2Env_NN_PEtHWS","2Env_NN_PEjHWS","2Env_NN_PEkHWS","2Env_NN_PEsHWS","2Env_NN_PEmHWS", "2Env_NN_PEnHWS"]
    #names = ["2Env_NN_PEAHWS","2Env_NN_PEEHWS","2Env_NN_PEVHWS","2Env_NN_PEHHWS","2Env_NN_PELHWS"]
    ## Extract data into something easily studied

    dims = [2,100]
    for name in names:
        for noise in [0]:#,0.25, 0.5,0.75,1,1.5,2,3,4,5]:
            for dim in dims:
                pass
                #main2_4(name,dim=dim,noise=noise,numE2=1)
                #main2_4_Plot(name,dim=dim,noise=noise,numE2=2)

            # Superpose 1 and 2 yamma
            #main2_4_PlotTogether(name,dim=[1,2],noise=noise,numE2=5)


    for name in names:
        pass
        #printAB(name)

    for cost in [0]:#,2,3,4]:
        for dim in dims:
            pass
            ## imshow with prop of cells in high. asymmetry vs frequency
            #figStrategiesRegion(names,cost,dim,noise=0) #imshow

            ## transition rate vs prop cells in high
            #figStrategiesPropTo(names,cost,dim=dim,noise=0) #lines

    for dim in [100]:
        pass
        ## asymmetry vs prop cells in high. Only for 2 and 0 env
        #figStrategiesAsymmetryo(names,'','log',dim=dim,noise=0)

    ## Difference in fitness sensing vs no sensing as a function of prop cells in high
    #figStrategiesDFitnessSensing(names,cost=[1,0],dims=2) #lognormal difference as propto in xcaxis

    for cost in [1]:
        pass
        ## Difference in fitness as a function of prop cells in high
        #figStrategiesEnvironmentEffect(names,noises = [0,0.25, 0.5,0.75,1,1.5,2,3,4,5],cost=cost,dims=[2,1]) #lognormal difference as propto in xcaxis

        ## Difference in fitness as a function fo the nois ein the environemnt. Only for two environemtns
        #figNoiseFitnessDifference(names,[0,1.,5],cost=cost,dims=[2,1])
        #figNoiseFitnessDifference(names,[0,0.25, 0.5,0.75,1,1.5,2,3,4,5],cost=cost,dims=[2,1])

        ## For 3 environments, difference in fitness as a function of the second environment. Only for optimum
        #figNoiseFitnessDifferenceOnlyOptimum3Env(names,cost=cost,dims=[2,1])

        # Noise only for optimum
        #figNoiseFitnessDifferenceOnlyOptimum([0,0.25, 0.5,0.75,1,1.5,2,3,4,5],cost=cost)

        #figGamma(names,dims=[2,1])


    ## For 3 environments, difference in fitness as a function of the second environment
    for cost in [0,1]:
        pass
        #fig3Strategies(names,cost=cost,dims=[2,1]) #lognor     mal difference as propto in xcaxis

    ## For 2/3 environments. Difference in fitness imshow
    for cost in [0,1]:
            pass
            ## imshow with prop of cells in high. asymmetry vs frequency
            #figStrategiesRegionDeltaFitness(names,cost,dims=[2,1])



  """

if __name__ == "__main__": main()
