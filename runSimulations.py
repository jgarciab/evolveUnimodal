import numpy as np
from numpy import random
import os
from scipy.stats import gamma, expon
import statsmodels.api as sm
import pylab as plt

class differential_evolution_optimizer(object):
  """
This is a python implementation of differential evolution
It assumes an evaluator class is passed in that has the following
functionality
data members:
 n              :: The number of parameters
 domain         :: a  list [(low,high)]*n
                   with approximate upper and lower limits for each parameter
 x              :: a place holder for a final solution

 also a function called 'target' is needed.
 This function should take a parameter vector as input and return a the function to be minimized.

 The code below was implemented on the basis of the following sources of information:
 1. http://www.icsi.berkeley.edu/~storn/code.html
 2. http://www.daimi.au.dk/~krink/fec05/articles/JV_ComparativeStudy_CEC04.pdf
 3. http://ocw.mit.edu/NR/rdonlyres/Sloan-School-of-Management/15-099Fall2003/A40397B9-E8FB-4B45-A41B-D1F69218901F/0/ses2_storn_price.pdf


 The developers of the differential evolution method have this advice:
 (taken from ref. 1)

If you are going to optimize your own objective function with DE, you may try the
following classical settings for the input file first: Choose method e.g. DE/rand/1/bin,
set the number of parents NP to 10 times the number of parameters, select weighting
factor F=0.8, and crossover constant CR=0.9. It has been found recently that selecting
F from the interval [0.5, 1.0] randomly for each generation or for each difference
vector, a technique called dither, improves convergence behaviour significantly,
especially for noisy objective functions. It has also been found that setting CR to a
low value, e.g. CR=0.2 helps optimizing separable functions since it fosters the search
along the coordinate axes. On the contrary this choice is not effective if parameter
dependence is encountered, something which is frequently occuring in real-world optimization
problems rather than artificial test functions. So for parameter dependence the choice of
CR=0.9 is more appropriate. Another interesting empirical finding is that rasing NP above,
say, 40 does not substantially improve the convergence, independent of the number of
parameters. It is worthwhile to experiment with these suggestions. Make sure that you
initialize your parameter vectors by exploiting their full numerical range, i.e. if a
parameter is allowed to exhibit values in the range [-100, 100] it's a good idea to pick
the initial values from this range instead of unnecessarily restricting diversity.

Keep in mind that different problems often require different settings for NP, F and CR
(have a look into the different papers to get a feeling for the settings). If you still
get misconvergence you might want to try a different method. We mostly use DE/rand/1/... or DE/best/1/... .
The crossover method is not so important although Ken Price claims that binomial is never
worse than exponential. In case of misconvergence also check your choice of objective
function. There might be a better one to describe your problem. Any knowledge that you
have about the problem should be worked into the objective function. A good objective
function can make all the difference.

Note: NP is called population size in the routine below.)
Note: [0.5,1.0] dither is the default behavior unless f is set to a value other then None.

  """

  def __init__(self,
               evaluator,
               population_size=50,
               f=None,
               cr=0.9,
               eps=1e-2,
               n_cross=1,
               max_iter=10000,
               monitor_cycle=200,
               out=None,
               show_progress=False,
               save_progress=False,
               show_progress_nth_cycle=1,
               insert_solution_vector=None,
               dither_constant=0.4,
                movAverageMutationRate = 0.,
                noise=0):

    self.movAverageMutationRate=movAverageMutationRate
    self.dither=dither_constant
    self.noise = noise
    self.show_progress=show_progress
    self.save_progress=save_progress
    self.show_progress_nth_cycle=show_progress_nth_cycle
    self.evaluator = evaluator
    self.population_size = population_size
    self.f = f
    self.cr = cr
    self.n_cross = n_cross
    self.max_iter = max_iter
    self.monitor_cycle = monitor_cycle
    self.vector_length = evaluator.n
    self.eps = eps
    self.population = []
    self.seeded = False
    if insert_solution_vector is not None:
      assert len( insert_solution_vector )==self.vector_length
      self.seeded = insert_solution_vector
    for ii in xrange(self.population_size):
      self.population.append( np.zeros(self.vector_length))


    self.scores = np.zeros(self.population_size) + 1000.
    self.optimize()
    self.best_score = np.min( self.scores )
    self.best_vector = self.population[( self.scores ).argmin() ]
    self.evaluator.x = self.best_vector


    if self.show_progress:
      self.evaluator.print_status(
            np.min(self.scores),
            np.mean(self.scores),
            self.population[ ( self.scores ).argmin() ],
            'Final')


  def optimize(self):
    # open file

    # initialise the population please
    self.make_random_population()
    # score the population please
    self.score_population()
    converged = False
    monitor_score = np.min( self.scores )
    self.count = 0
    cx = 0
    while not converged:
      self.evolve()
      location = (self.scores).argmin()
      if self.show_progress:
        if self.count%self.show_progress_nth_cycle==0:
          # make here a call to a custom print_status function in the evaluator function
          # the function signature should be (min_target, mean_target, best vector)
          self.evaluator.print_status(
            np.min(self.scores),
            np.mean(self.scores),
            self.population[ ( self.scores ).argmin() ],
            self.count)
      if self.save_progress:
        self.evaluator.fname.write("%d, %f, %f" %(self.count,np.min(self.scores),np.mean(self.scores)))
        for item in self.population[ ( self.scores ).argmin() ]:
          self.evaluator.fname.write(", %e" % item)
        if self.count%20==0:
            print self.count, self.evaluator.fname.name, np.min(self.scores), self.population[ ( self.scores ).argmin() ]
            #print self.count
            #vector = self.population[ ( self.scores ).argmin()][:-1]
            #x = np.linspace(0.01, 100., num=100) # values for x-axis
            #d = np.zeros(100)

            #for jj in range(0,len(vector)-1,3):
                #d += vector[jj]*gamma.pdf(x, vector[jj+1], loc=0, scale=vector[jj+2]) # probability distribution
            #plt.plot(d)
            #plt.show()


        self.evaluator.fname.write("\n")

      self.count += 1
      if self.count%self.monitor_cycle==0:
        if (monitor_score-np.min(self.scores) ) < self.eps:
          converged = True
        else:
         monitor_score = np.min(self.scores)
      rd = (np.mean(self.scores) - np.min(self.scores) )
      rd = rd*rd/(np.min(self.scores)*np.min(self.scores) + self.eps )


      if ( rd < self.eps):
        cx += 1

      if self.count>=self.max_iter :
        converged = True

      if cx > 20:
          converged = True


    if self.save_progress:
        self.evaluator.fname.close()

    return None

  def make_random_population(self):
    for ii in xrange(self.vector_length):
      delta  = self.evaluator.domain[ii][1]-self.evaluator.domain[ii][0]
      offset = self.evaluator.domain[ii][0]
      random_values = np.random.random(self.population_size)
      random_values = random_values*delta+offset
      # now please place these values ni the proper places in the
      # vectors of the population we generated
      for vector, item in zip(self.population,random_values):
        vector[ii] = item
    if self.seeded is not False:
      self.population[0] = self.seeded

    self.upper_bound = np.asarray([_[1] for _ in self.evaluator.bounder])
    self.lower_bound = np.asarray([_[0] for _ in self.evaluator.bounder])

    """
    for vector in self.population:
        x = np.linspace(0.01, 100., num=100) # values for x-axis
        d = np.zeros(100)
        for jj in range(0,len(vector)-1,3):
            d += vector[jj]*gamma.pdf(x, vector[jj+1], loc=0, scale=vector[jj+2]) # probability distribution
        d /= np.sum(d)
        plt.plot(d)

    plt.show()
    """

  def score_population(self):
    for ii,vector in enumerate(self.population):
      tmp_score = self.evaluator.target(vector,0)
      self.scores[ii]=tmp_score

  def evolve(self):
    #print self.scores[(self.scores ).argmin()]
    for ii in xrange(self.population_size):
      if self.noise != 0:
        self.scores[ii] = self.evaluator.target( self.population[ii],self.count )
      np.random.seed()
      permut = np.random.permutation(self.population_size-1)
      # make parent indices
      i1=permut[0]
      if (i1>=ii):
        i1+=1
      i2=permut[1]
      if (i2>=ii):
        i2+=1
      i3=permut[2]
      if (i3>=ii):
        i3+=1
      """
      x1 = self.population[ i1 ]
      x2 = self.population[ i2 ]
      x3 = self.population[ i3 ]

      if self.f is None:
        use_f = random.random()/2.0 + 0.5
      else:
        use_f = self.f

      vi = x1 + use_f*(x2-x3)
      # crossover
      mask = np.random.random(self.vector_length)
      test_vector = (mask < 0.9)*vi + (mask>0.9)*self.population[ii]
      test_vector[test_vector<self.lower_bound] = self.lower_bound[test_vector<self.lower_bound]
      test_vector[test_vector>self.upper_bound] = self.upper_bound[test_vector>self.upper_bound]
      """

      if self.count < 50 or np.random.random()<0.8:
        x1 = self.population[ i1 ]#self.population[ i1 ]#
      else:
        x1 = self.population[ ( self.scores ).argmin()]#self.population[ i1 ]#self.population[ i1 ]#
      x2 = self.population[ i2 ]
      x3 = self.population[ i3 ]

      if self.f is None:
        use_f = random.random()/2.0 + 0.5
      else:
        use_f = self.f

      vi = x1 + use_f*(x2-x3)
      # crossover
      mask = np.random.random(self.vector_length)
      test_vector = (mask < 0.9)*vi + (mask>0.9)*self.population[ii]
      test_vector[test_vector<self.lower_bound] = self.lower_bound[test_vector<self.lower_bound]
      test_vector[test_vector>self.upper_bound] = self.upper_bound[test_vector>self.upper_bound]

      # moving average
      if np.random.random() < self.movAverageMutationRate:
          rN = 3#np.random.randint(2,5)*2-1
          t1,t2= np.sum(test_vector[:40]),np.sum(test_vector[40:-1])
          test_vector = np.concatenate([test_vector[:rN/2], (np.convolve(test_vector[:-1]**rN, np.ones((rN,))/float(rN), mode='valid'))**rN,test_vector[(-rN-1)/2:-1]**rN,[test_vector[-1]]])
          test_vector[:40] /= np.sum(test_vector[:40]) / t1
          test_vector[40:-1] /= np.sum(test_vector[40:-1]) / t2

      if np.random.random() < self.movAverageMutationRate:
          if random.random() < 0.5:
            test_vector[:40] = 1./2 * (test_vector[:40]+ test_vector[1:41])
            test_vector[40:-2] = 1./2 * (test_vector[41:-1]+ test_vector[40:-2])

          else:
            test_vector[:40] = 1./2 * (test_vector[:40]+ test_vector[1:41])
            test_vector[41:-1] = 1./2 * (test_vector[41:-1]+ test_vector[40:-2])

      if np.random.random() < self.movAverageMutationRate:
          if random.random() < 0.5:
            test_vector[:40] *= 1.01
          else:
            test_vector[40:-1] *= 1.01

      # bounder
      test_score = self.evaluator.target( test_vector,self.count )

      if test_score < self.scores[ii]:
        self.scores[ii] = test_score
        self.population[ii] = test_vector



  def show_population(self):
    print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    for vec in self.population:
      print list(vec)
    print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"


class Init(object):
    def __init__(self,evaluator, suddenness,numChanges,args,dim,noise = 0):
        evaluator.numEnv = int(args[0])


        if noise == 0 or evaluator.numEnv == 0:
            y = [0]*1000
            if evaluator.numEnv == 0:
                lenIter = 1000
            else:
                lenIter = 2
        else:
            x = np.linspace(0.0, 100., num=101)
            tt  =expon.pdf(x,scale=noise,loc=0)
            tt = tt/np.sum(tt)
            if evaluator.numEnv == 2:
                lenIter = 200
            else:
                lenIter = 50

            y = []
            for i,t in enumerate(tt):
                y += [int(x[i])]*int(lenIter*2*t)

        evaluator.noise = y

        costArr = ['0','0','0.01','0.03','0.1']
        cost = costArr[suddenness]

        if evaluator.numEnv == 0:
            a = float(args[args.find("0Env_")+5] + "." + args[args.find("0Env_")+6:-2])
            j = 1.5/a

            np.random.seed(2+0)
            x = np.linspace(0.01, 100., num=101)
            tt  =gamma.pdf(x,a,scale=j,loc=0)
            tt = tt/np.sum(tt)
            y = []
            for i,t in enumerate(tt):
                y += [int(11*x[i])]*int(1000*t)

            evaluator.env = np.random.choice([int(_) for _ in y],size=len(y),replace=False)
            print set(evaluator.env)

            evaluator.trajectory = dict()
            i = 0
            for s in range(len(evaluator.env)):
                i += int(10000/numChanges)
                evaluator.trajectory[i] = s

        if evaluator.numEnv == 1:
            s = int(args[args.find("0Env_")+6:-2])
            print(s)

            evaluator.env = [s,s]
            if 1:
                evaluator.trajectory = dict()
                evaluator.trajectory[1000] = 0



        elif evaluator.numEnv == 2:
            evaluator.env = [0,100]

            if args[-4] == 'A': x2 = 0.999999 #1000000
            elif args[-4] == 'B': x2 = 0.999998 #1000000
            elif args[-4] == 'C': x2 = 0.999995 #1000000

            elif args[-4] == 'E': x2 = 0.99999 #100000
            elif args[-4] == 'F': x2 = 0.99998 #50000
            elif args[-4] == 'G': x2 = 0.99995 #20000

            elif args[-4] == 'V': x2 = 0.9999 #10000
            elif args[-4] == 'W': x2 = 0.9998 #5000
            elif args[-4] == 'X': x2 = 0.9995 #2000

            elif args[-4] == 'H': x2 = 0.999 #1000

            elif args[-4] == 'I': x2 = 0.9960#80 #500
            elif args[-4] == 't': x2 = 0.9958#79 #400
            elif args[-4] == 'j': x2 = 0.9956#78 #333
            elif args[-4] == 'k': x2 = 0.9954#77 #434
            elif args[-4] == 's': x2 = 0.9952#76 #434
            elif args[-4] == 'm': x2 = 0.9950#75 #434
            elif args[-4] == 'n': x2 = 0.9948#74 #434


            #elif args[-4] == 'I': x2 = 0.9980#56#80 #500
            #elif args[-4] == 't': x2 = 0.9979#54#79 #400
            ##elif args[-4] == 'j': x2 = 0.9978#52#78 #333
            #elif args[-4] == 'k': x2 = 0.9977#50#77 #434
            #elif args[-4] == 's': x2 = 0.9976#48#76 #434
            #elif args[-4] == 'm': x2 = 0.9975#46#75 #434
            #elif args[-4] == 'n': x2 = 0.9974#44#74 #434


            elif args[-4] == 'o': x2 = 0.9973 #434
            elif args[-4] == 'p': x2 = 0.9972 #434
            elif args[-4] == 'q': x2 = 0.9971 #434
            elif args[-4] == 'r': x2 = 0.997 #434
            elif args[-4] == 'J': x2 = 0.995 #200



            elif args[-4] == 'L': x2 = 0.99 #100

            if args[-3] == 'V': x3 = 0.9999
            elif args[-3] == 'H': x3 = 0.999
            elif args[-3] == 'L': x3 = 0.99
            elif args[-3] == 'A': x3 = 0.999999 #1000000

            if args[-6] == 'P':
                evaluator.trajectory = dict()
                s = 1
                i = 0
                while(len(evaluator.trajectory)<lenIter):
                    if s == 0:
                        #v5 (Very low freq in High stress)
                        i += int(np.ceil(1000.*1./(1-x2)/numChanges))
                    else:
                        i += int(np.ceil(1000.*1./(1-x3)/numChanges))
                    s = (s-1)*(-1)
                    evaluator.trajectory[i] = s


        elif evaluator.numEnv == 3:
            evaluator.env = [0,11,100]

            if args[-5] == 'A': x1 = 0.999999 #1000000
            elif args[-5] == 'E': x1 = 0.99999 #100000
            elif args[-5] == 'V': x1 = 0.9999 #10000
            elif args[-5] == 'H': x1 = 0.999 #1000
            elif args[-5] == 'L': x1 = 0.99 #100

            if args[-4] == 'A': x2 = 0.999999 #1000000
            elif args[-4] == 'E': x2 = 0.99999 #100000
            elif args[-4] == 'V': x2 = 0.9999 #10000
            elif args[-4] == 'H': x2 = 0.999 #1000
            elif args[-4] == 'L': x2 = 0.99 #100

            if args[-3] == 'H': x3 = 0.999

            if args[-7] == 'P':
                #Regular
                evaluator.trajectory = dict()
                envOrder = [0,1,0,2]

                s = 1
                i = 0
                while(len(evaluator.trajectory)<2*lenIter):
                    if envOrder[s%4] == 1:
                        i += int(np.ceil(1./(1-x2)/numChanges))
                    elif envOrder[s%4] == 2:
                        i += int(np.ceil(1./(1-x3)/numChanges))
                    else:
                        i += int(0.5*np.ceil(1./(1-x1)/numChanges))

                    s+=1
                    evaluator.trajectory[i] = envOrder[s%4]



        if args[-2] == 'S':
            evaluator.arrayCost = []
            evaluator.arrayCost.append(np.loadtxt('allCostsSt_S'+'0'+'.txt'))
            evaluator.arrayCost.append(np.loadtxt('allCostsSt_S'+cost+'.txt'))
            evaluator.selection = 1

        elif args[-2] == 'W':
            evaluator.arrayCost = []
            evaluator.arrayCost.append(np.loadtxt('allCostsSt_W'+'0'+'.txt'))
            evaluator.arrayCost.append(np.loadtxt('allCostsSt_W'+cost+'.txt'))
            evaluator.selection = 0

        else:
            print "Finish with SS or WS"
            raise
        evaluator.optVal = [evaluator.arrayCost[1][:,i].argmax() for i in range(101)]
        evaluator.gamma1Env = np.loadtxt("gamma1EnvOptimum.txt")

        ## Global variables
        evaluator.sud = suddenness
        evaluator.trajectoryX = evaluator.trajectory
        evaluator.trajectory = sorted([_ for _ in evaluator.trajectory])

        print evaluator.trajectoryX


class EvolveNoiseFromHistLogNormal(object):
    def __init__(self, suddenness,numChanges,args,dim,noise = 0):
        self.fname = open("./dataDE/"+str(noise)+args+str(dim)+str(suddenness)+str(numChanges)+"2obs.txt","w")
        Init(self,suddenness,numChanges,args,dim,noise)

        self.x = None
        self.n = dim*3+1
        self.dim = dim*3+1
        if dim == 1:
            self.domain = [(0.,1.), (0.5,100.),(10.,400.)] + [(0,1)]
            self.bounder = [(0.,10.), (0.5,100),(10.,4000.)] +[(0,1)]
        else:
            if dim %2 != 0:
                 dim -= 1
                 print "Dimensions reduced"
            self.domain = [(0.,1.), (0.5,2),(10.,400.),(0.,1.), (2,100),(10.,400.)]*(dim/2) + [(0,1)]
            self.bounder = [(0.,10.), (0.5,100),(10,4000.),(0.,10.), (0.5,100),(10.,4000.)]*(dim/2) + [(0,1)]

        self.optimizer =  differential_evolution_optimizer(self,max_iter=500 ,population_size=400,
                                                           n_cross=1,cr=0.9, eps=1e-15, show_progress=False,save_progress=True,noise=noise)


    def target(self, vector,seed):
        random.seed(100*seed+0)
        x = np.linspace(0.01, 10000., num=100) # values for x-axis
        d = np.zeros(100)

        w = 0
        for jj in range(0,len(vector)-1,3):
            d += vector[jj]*gamma.cdf(x, vector[jj+1], loc=0, scale=vector[jj+2]) # probability distribution
            w += vector[jj]

        d = np.diff(np.concatenate([[0],d]))
        sense = np.round(vector[-1])
        timePointAll = d/w


        timePoint = np.copy(timePointAll)

        currEnv = 1
        sumT = 0
        prevchange = 0

        np.random.shuffle(self.noise)
        for i,change in enumerate(self.trajectory):
            if currEnv == 0:
                env = self.env[currEnv] + self.noise[i]
                temp = np.copy(timePointAll)
            else:
                env =  self.env[currEnv] - self.noise[i]
                a,b = self.gamma1Env[:,env]
                temp = np.diff(np.concatenate([[0],gamma.cdf(x, a, loc=0, scale=b)]))# probability distribution


            if sense == 1:
                opt = self.arrayCost[1][:,env]
            else:
                opt = self.arrayCost[0][:,env]

            inter = change-prevchange
            #print "1",i,currEnv,env,inter,change

            prevchange = change
            if sense == 0 or self.sud == 0:
                growth = np.sum(timePoint[opt>-1]*2**opt[opt>-1])
                if growth == 0: return 1.
                sumT +=  1.*inter*np.log2(growth)

            else:
                t2 = temp
                #First see who grows
                growth = np.sum(timePoint[opt>-1]*2**opt[opt>-1])
                if growth == 0: return 1.

                #Now switch. Fast changes
                sumT += 1.*np.log2(growth)
                sumT += 1.*(inter-1)*np.log2(np.sum(t2[opt>-1]*2**opt[opt>-1]))
                #print 1.*np.log(growth),1.*(inter-1)*np.log(np.sum(t2 + t2 * opt))

            currEnv = self.trajectoryX[change]
            #print "2",i,currEnv,env,inter,change

        fitness = sumT/self.trajectory[-1]#np.exp(sumT/self.trajectory[-1])-1.
        #print fitness

        if 0:
            penalty = 0.1*np.sum(np.abs(np.diff(timePointAll))>0.01) #0.1 for each sudden change in concentration
            fitness =  fitness-penalty
        else:
            fitness =  fitness

        if np.isnan(fitness): return 2.
        else: return -fitness


    def print_status(self, mins,means,vector,txt):
        print txt,mins, means, list(vector)

class EvolveNoiseFromHistStd(object):
    def __init__(self, suddenness,numChanges,args,dim,noise = 0):
        Init(self,suddenness,numChanges,args,dim,noise)
        self.fname = open("./dataDE/"+str(noise)+args+str(dim)+str(suddenness)+str(numChanges)+"2STDobs.txt","w")


        self.x = None
        self.n = 101
        self.dim = 101

        self.domain = [(0.,1.)] *100 + [(0,1)]
        self.bounder = [(0.,1.)] *100 + [(0,1)]

        self.optimizer =  differential_evolution_optimizer(self,max_iter=500 ,population_size=500,
                                                           n_cross=1,cr=0.9, eps=1e-15, show_progress=False,
                                                           save_progress=True,movAverageMutationRate = 0.1 ,noise=noise)


    def target(self, vector,seed):
        random.seed(100*seed+0)
        d = vector[:-1]

        sense = np.round(vector[-1])
        timePointAll = d/np.sum(d)
        timePoint = np.copy(timePointAll)

        currEnv = 1
        sumT = 0
        prevchange = 0


        np.random.shuffle(self.noise)
        for i,change in enumerate(self.trajectory):
            if currEnv == 0:
                env = self.env[currEnv] + self.noise[i]
                temp = np.copy(timePointAll)
            else:
                env =  self.env[currEnv] - self.noise[i]
                temp = np.zeros(100)
                temp[self.optVal[env]] = 1.


            if sense == 1:
                opt = self.arrayCost[1][:,env]
            else:
                opt = self.arrayCost[0][:,env]

            inter = change-prevchange
            #print inter, envX[currEnv]
            prevchange = change

            if sense == 0 or self.sud == 0:
                growth = np.sum(timePoint[opt>-1]*2**opt[opt>-1])
                if growth == 0: return 1.
                sumT +=  1.*inter*np.log2(growth)

            else:
                t2 = temp
                #First see who grows
                growth = np.sum(timePoint[opt>-1]*2**opt[opt>-1])
                if growth == 0: return 1.

                #Now switch. Fast changes
                sumT += 1.*np.log2(growth)
                sumT += 1.*(inter-1)*np.log2(np.sum(t2[opt>-1]*2**opt[opt>-1]))
                #print 1.*np.log(growth),1.*(inter-1)*np.log(np.sum(t2 + t2 * opt))


            currEnv = self.trajectoryX[change]


        #fitness = np.exp(sumT/self.trajectory[-1])-1.
        fitness = sumT/self.trajectory[-1]

        if 0:
            penalty = 0.1*np.sum(np.abs(np.diff(timePointAll))>0.01) #0.1 for each sudden change in concentration
            fitness =  fitness-penalty
        else:
            fitness =  fitness

        if np.isnan(fitness): return 2.
        else: return -fitness


    def print_status(self, mins,means,vector,txt):
        print txt,mins, means, list(vector)



def run(pF):
    import time
    random.seed(64+0)

    if pF[3] == 100:
        fname = str(pF[4])+pF[2]+str(pF[3])+str(pF[0])+str(pF[1])+"2STDobs.txt"
    else:
        fname = str(pF[4])+pF[2]+str(pF[3])+str(pF[0])+str(pF[1])+"0obs.txt"


    if fname in os.listdir('./dataDE/'):
        print fname, os.path.getsize('./dataDE/'+fname)
        if os.path.getsize('./dataDE/'+fname) > 1000:
            print time.ctime(os.path.getmtime('./dataDE/'+fname))
            pass#return None

    if pF[3] == 100:
        EvolveNoiseFromHistStd(pF[0],pF[1],pF[2],dim=pF[3],noise=pF[4])
    else:
        EvolveNoiseFromHistLogNormal(pF[0],pF[1],pF[2],dim=pF[3],noise=pF[4])
    #




def main():
  from multiprocessing import Pool


    #de2((sudden,numChanges,tit,0))
  #names =   ["0Env_00010SS","0Env_000250SS","0Env_00050SS","0Env_00750SS","0Env_00250SS","0Env_0050SS","0Env_00750SS"]#,"0Env_030SS","0Env_050SS","0Env_070SS","0Env_090SS","0Env_100SS","0Env_020SS","0Env_040SS","0Env_060SS","0Env_080SS"]



           #"2Env_NN_PEIHWS","2Env_NN_PEtHWS","2Env_NN_PEjHWS","2Env_NN_PEkHWS","2Env_NN_PEsHWS","2Env_NN_PEmHWS",
           #"2Env_NN_PEnHWS", "2Env_NN_PEoHWS","2Env_NN_PEpHWS","2Env_NN_PEqHWS","2Env_NN_PErHWS"]

  #names = ["2Env_NN_PEIHWS","2Env_NN_PEtHWS","2Env_NN_PEjHWS","2Env_NN_PEkHWS","2Env_NN_PEsHWS","2Env_NN_PEmHWS", "2Env_NN_PEnHWS"]






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
  #names = ["2Env_NN_PEAHSS","2Env_NN_PEEHSS","2Env_NN_PEVHSS","2Env_NN_PEHHSS","2Env_NN_PELHSS"]


  #names = ["2Env_NN_PEAHWS","2Env_NN_PEEHWS","2Env_NN_PEVHWS","2Env_NN_PEHHWS","2Env_NN_PELHWS"]
  names = ["2Env_NN_PEIHWS","2Env_NN_PEtHWS","2Env_NN_PEjHWS","2Env_NN_PEkHWS","2Env_NN_PEsHWS","2Env_NN_PEmHWS", "2Env_NN_PEnHWS"]
  possibleFactors = []
  for numChanges in  [1,3,10,30,100]:
      for sudden in range(1):
          for dim in [100]:
              for noise in [0]:
                  for name in names:
                    #pF = (sudden,numChanges,name,dim,noise)
                    #EvolveNoiseFromHistLogNormal(pF[0],pF[1],pF[2],dim=pF[3],noise=pF[4])
                    possibleFactors.append((sudden,numChanges,name,dim,noise))


  """
  possibleFactors = []
  for dim in [100]:
    for noise in [0]:
        for end in ["A","E","V","H","L"]:
            possibleFactors.append((0,10,"3Env_0102_PEA"+end+"HSS",dim,noise))
            possibleFactors.append((0,10,"3Env_0102_PEE"+end+"HSS",dim,noise))
            #possibleFactors.append((0,10,"3Env_0102_PEV"+end+"HSS",dim,noise))
            #possibleFactors.append((0,10,"3Env_0102_PEH"+end+"HSS",dim,noise))
            #possibleFactors.append((0,10,"3Env_0102_PEL"+end+"HSS",dim,noise))

            #possibleFactors.append((1,100,"3Env_0102_PEA"+end+"HSS",dim,noise))
            #possibleFactors.append((1,100,"3Env_0102_PEE"+end+"HSS",dim,noise))
            #possibleFactors.append((1,100,"3Env_0102_PEV"+end+"HSS",dim,noise))
            #possibleFactors.append((1,10,"3Env_0102_PEH"+end+"HSS",dim,noise))
            #possibleFactors.append((1,1,"3Env_0102_PEL"+end+"HSS",dim,noise))

  """



  """

  names = ["2Env_NN_PEAHSS","2Env_NN_PEEHSS","2Env_NN_PEVHSS","2Env_NN_PEHHSS","2Env_NN_PELHSS"]
  possibleFactors = []
  for dim in [1,2]:
    for noise in [0.25, 0.5,0.75,1,1.5,2,3,4,5]:
        #pF = (1,100,"2Env_NN_PEAHSS",dim,noise)
        #EvolveNoiseFromHistLogNormal(pF[0],pF[1],pF[2],dim=pF[3],noise=pF[4])

        # Not run for 1, run for 0
        for name in names:
            possibleFactors.append((0,10,name,dim,noise))

        #Run already for 0 and 1
        possibleFactors.append((1,100,"2Env_NN_PEAHSS",dim,noise))
        possibleFactors.append((1,100,"2Env_NN_PEEHSS",dim,noise))
        possibleFactors.append((1,100,"2Env_NN_PEVHSS",dim,noise))
        possibleFactors.append((1,10,"2Env_NN_PEHHSS",dim,noise))
        possibleFactors.append((1,1,"2Env_NN_PELHSS",dim,noise))

  """

  """
  for stress in range(100,101):
      if stress < 10:
          s = "0"+str(stress)
      else:
          s = str(stress)
      name = "1Env_"+s+"SS"
      pF =(0,1,name,1,0)
      EvolveNoiseFromHistLogNormal(pF[0],pF[1],pF[2],dim=pF[3],noise=pF[4])
  """
  """
  names = ["2Env_NN_PEIHWS","2Env_NN_PEtHWS","2Env_NN_PEjHWS","2Env_NN_PEkHWS","2Env_NN_PEsHWS","2Env_NN_PEmHWS", "2Env_NN_PEnHWS"]
  possibleFactors = []
  for numChanges in  [1,3,10,30,100]:
      for sudden in range(1):
          for dim in [100]:
              for noise in [0]:
                  for name in names:
                    #pF = (sudden,numChanges,name,dim,noise)
                    #EvolveNoiseFromHistLogNormal(pF[0],pF[1],pF[2],dim=pF[3],noise=pF[4])
                    possibleFactors.append((sudden,numChanges,name,dim,noise))

  """
  pool = Pool(processes=8)
  pool.map(run, possibleFactors)
  pool.close()
  pool.join() #zombie processes without this, will fill up memory
  print "OK"


if __name__ == "__main__":
  main()
  #EvolveNoiseFromHistStd(1,1,"2Env_NN_PEVHSS",dim=100,noise=0)
