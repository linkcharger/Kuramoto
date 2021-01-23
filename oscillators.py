# %% setup & construction %%
import math
import random as rand
import numpy as np
import matplotlib.pyplot as plt
from numba import jit, double
# from concurrent.futures import ThreadPoolExecutor
from copy import deepcopy



############################################### options for optimisation #############################################
#
# jit
# 
# double sum -> even better, average phase, then ...? wait what?
# 
#####################################################################################################################





class Oscillator:

    def __init__(self, initPhase, omega):
        self.omega = omega
        self.lastTheta = 0.0                                            # gets value immediately, but instantiate type
        self.currentTheta = initPhase                                   # handed over immediately on first run to lastTheta




class OscPopulation:

    def __init__(self, omegaDistr = 'normal'):
        # list of objects  
        self.list_os = np.zeros(N, dtype = Oscillator)                      # pre-allocate with fixed size -> speed gain?
        self.omegaDistr = omegaDistr
       
       # complete initialisation
        for n in range(N):
            if omegaDistr == 'normal': omega = rand.normalvariate(0,1)
            elif omegaDistr == 'uniform': omega = rand.uniform(-0.5, 0.5)
                    
            initPhase = rand.uniform(0, math.tau)        
            
            self.list_os[n] = Oscillator(initPhase, omega)

        print('oscillators initialised with', omegaDistr.upper(), 'distibution of omegas' )



    def reRandomise(self):

        for n in range(N):
            # omega
            if self.omegaDistr == 'normal': newOmega = rand.normalvariate(0.0,1.0)
            elif self.omegaDistr == 'uniform': newOmega = rand.uniform(-0.5, 0.5)  
            # theta
            newInitPhase = rand.uniform(0, math.tau)      
            
            self.list_os[n].omega = newOmega
            self.list_os[n].lastTheta = 0.0                                             # gets value immediately, but instantiate type
            self.list_os[n].currenTtheta = newInitPhase                                    # handed over immediately on first run to lastTheta
            

        print('3 variables re-randomised' )




    # def changeThetas(self):

    #     for n in range(N):
    #         initPhase = rand.uniform(0, math.tau)

    #         # update thetas
    #         self.list_os[n].lastTheta = 0.0
    #         self.list_os[n].currentTheta = initPhase
    #         # dont touch omegas


    # def changeOmegas(self):

    #     for n in range(N):
    #         omega = rand.uniform(-0.5, 0.5)

    #         # update omegas
    #         self.list_os[n].omega = omega
    #         # dont touch thetas

        



    ############ calculate r ############
    def calc_r(self):                                       
        
        ###### one way to do it #######
        # complex sum of all phases 
        # realSum = 0 
        # imaginarySum = 0

        # for n in range(N): 
        #     realSum += math.cos(self.list_os[n].theta[-1])
        #     imaginarySum += math.sin(self.list_os[n].theta[-1])
        
        # realSum = realSum/N 
        # imaginarySum = imaginarySum/N
        # r = math.sqrt(realSum * realSum + imaginarySum * imaginarySum)


        ###### another way to do it #######
        sum = 0.0 + 0.0j
        for n in range(N):
            sum += math.e**(1.0j * self.list_os[n].currentTheta)                             
        sum = sum / N 

        # r = math.sqrt(sum.real**2 + sum.imag**2)
        r = abs(sum)
        
        return r

    

    # wrapper method
    def stepAll(self, K):

        # getter: from objects
        lastThetaList = []
        currentThetaList = []
        omegaList = []

        for n in range(N):
            lastThetaList.append(self.list_os[n].lastTheta)
            currentThetaList.append(self.list_os[n].currentTheta)
            omegaList.append(self.list_os[n].omega)

        # converting
        lastThetaArray = np.array(lastThetaList)
        currentThetaArray = np.array(currentThetaList)
        omegaArray = np.array(omegaList)

        # up until here its totally fine, values are correctly transmitted

        # calculating, handing over
        newLastThetaArray, newCurrentThetaArray = self._stepAll(lastThetaArray, currentThetaArray, omegaArray, K)

        # conversion, just in case
        lastThetaList = newLastThetaArray.tolist()
        currentThetaList = newCurrentThetaArray.tolist()



        # setter: back to OOP
        for n in range(N):
            self.list_os[n].lastTheta = lastThetaList[n]
            self.list_os[n].currentTheta = currentThetaList[n]

    
    
    # step FOR time t -> use t-1 ####################################
    @staticmethod
    @jit(nopython = True)            #('double[:](double[:], double[:], double[:], int8)', nopython = True)
    def _stepAll(lastThetaArray, currentThetaArray, omegaArray, K):
        
        for n in range(N):                                                                  # step through time -> hand over value
            lastThetaArray[n] = currentThetaArray[n]

        # calculate new thetas
        for n in range(N):
            
            sum = 0.0
            for j in range(N):                                                              # calculate differential sum for ONE oscillator n for time t
                sum += math.sin(lastThetaArray[j] - lastThetaArray[n])                  

            theta_dot_t = omegaArray[n] + K/N * sum                                 # theta_dot_t for oscillator n
            
            theta_t = lastThetaArray[n] + dt * theta_dot_t                          # new theta for oscillator n -> euler step

            currentThetaArray[n] = theta_t                                          # going down list of objects, pick object, dial into theta list, append

        return lastThetaArray, currentThetaArray




    def runK(self):
        
        r_critList= []
                        
        for K in K_range:  
            self.reRandomise()

            for t in t_range:
                self.stepAll(K)

            r_crit = self.calc_r()
            r_critList.append(r_crit)
            print('K = ' + str(round(K,3)) + ' | r_crit = ' + str(round(r_crit, 3)))
            
        return r_critList



    def runT(self):
        if len(K_range) > 1:
            rList = [[]]

            for K in K_range:  
                print('\nK =', K)
                self.reRandomise()
                # print('(reRandomised called by K loop)')

                rList.append([])

                for t in t_range:
                    # step all oscillators forward INTO timeperiod t
                    self.stepAll(K)

                    # calculate the population's coherence at each time step
                    r = self.calc_r()
                    rList[K].append(r)
                            
                    # status
                    if t % 10 == 0: print('t =', t, 'done') 

        elif len(K_range) == 1:
            rList = []

            for K in K_range:  
            
                for t in t_range:
                    # step all oscillators forward INTO timeperiod t
                    self.stepAll(K)

                    # calculate the population's coherence at each time step
                    r = self.calc_r()
                    rList.append(r)
                            
                    # status
                    if t % 10 == 0: print('t =', t, 'done') 
        
        
        return rList
                







# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% task 1: normal omegas, K-vs-r %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%time


N = 200             # 1000  -  100 took about 45min
T = 100             # 100
dt = 0.01
K = 4
dk = 0.1

numberOfTimes = int(T/dt)
numberOfK = int(K/dk)
t_range = [round(i * dt, 4) for i in range(numberOfTimes + 1)]
K_range =  [round(i * dk, 4) for i in range(numberOfK + 1)]            


pop1 = OscPopulation()
r_critList = pop1.runK() 


# graphics 
plt.figure(figsize = (10,6))
plt.title('System state diagram')
plt.xlim(0, 4)
plt.ylim(0, 1)
plt.xlabel('coupling strength K')
plt.ylabel('coherence r')

plt.plot(K_range, r_critList, 'ko')
filename = 'graphics/1 K-vs-r_' + 'omegaDistr=' + pop1.omegaDistr + '_N=' + str(N) + '.pdf'
plt.savefig(filename, dpi = 200, bbox_inches = 'tight')
plt.show()








# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% task 2: normal omegas, t-vs-r %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%time


N = 500             # 1000
T = 100             # 100
dt = 0.01

numberOfTimes = int(T/dt)
t_range = [round(i * dt, 4) for i in range(numberOfTimes + 1)]      
K_range = [1, 2, 3]

pop2 = OscPopulation('normal')
rList = pop2.runT()

### graphing ###
plt.figure(figsize = (10,6))
plt.title('Evolution of r')
plt.ylim(0, 1)
plt.xlabel('time t')
plt.ylabel('coherence r')

for K in K_range:
    plt.plot(t_range, rList[K], label = 'K = ' + str(K))
plt.legend()
filename = 'graphics/2 t-vs-r_' + 'omegaDistr=' + pop2.omegaDistr + '_N=' + str(N) + '.pdf'
plt.savefig(filename, dpi = 200, bbox_inches = 'tight')
plt.show()








# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% task 3: uniform omegas, K-vs-r %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%time


N = 400                   # 2000    - 100, 100, took 6min
T = 200                   # 200
dt = 0.05
K = 1.5
dk = 0.03

numberOfTimes = int(T/dt)
numberOfK = int(K/dk)
t_range = [round(i * dt, 4) for i in range(numberOfTimes + 1)]
K_range =  [round(i * dk, 4) for i in range(numberOfK + 1)]  

pop3 = OscPopulation('uniform')
r_critList = pop3.runK() 


# graphics 
plt.figure(figsize = (10,6))
plt.title('System state diagram')
plt.xlim(0, 4)
plt.ylim(0, 1)
plt.xlabel('coupling strength K')
plt.ylabel('coherence r')

plt.plot(K_range, r_critList, 'ko')
filename = 'graphics/3 K-vs-r_' + 'omegaDistr=' + pop3.omegaDistr + '_N=' + str(N) + '.pdf'
plt.savefig(filename, dpi = 200, bbox_inches = 'tight')
plt.show()


# graphics 
plt.figure(figsize = (10,6))
plt.title('System state diagram')
plt.ylim(0, 1)
plt.xlabel('coupling strength K')
plt.ylabel('coherence r')

plt.plot(K_range, r_critList, 'ko')
filename = 'graphics/3 K-vs-r_' + 'omegaDistr=' + pop3.omegaDistr + '_zoomed' + '_N=' + str(N) + '.pdf'
plt.savefig(filename, dpi = 200, bbox_inches = 'tight')
plt.show()





# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% task 4: uniform omegas, fixed omegas, change thetas, t-vs-r repeated %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%time
# --> same natural frequencies as first time, but start at different positions every time

N = 400                      # 2000
T = 200                      # 200
dt = 0.05
runs = 10 

numberOfTimes = int(T/dt)
t_range = [round(i * dt, 4) for i in range(numberOfTimes + 1)]
K_range = [2]                # 1

fixedOmegaPop = OscPopulation('uniform')
# savedOmegas = deepcopy(fixedOmegaPop)
rLists = []


for run in range(runs):
    _fixedOmegaPop = OscPopulation('uniform')

    for n in range(N): 
        _fixedOmegaPop.list_os[n].omega = fixedOmegaPop.list_os[n].omega
    print('original omegas restored')

    rNewList = _fixedOmegaPop.runT()
    rLists.append(rNewList)

    print('############# run', str(run), 'done #############\n') 


### graphing ###
plt.figure(figsize = (10,6))
plt.title('Evolution of r')
plt.ylim(0,1)
plt.xlabel('time t')
plt.ylabel('coherence r')

for run in range(runs):
    plt.plot(t_range, rLists[run])

filename = 'graphics/4 t-vs-r' + '_fixedOmegas' + '_omegaDistr=' + fixedOmegaPop.omegaDistr + '_N=' + str(N) + '.pdf'
plt.savefig(filename, dpi = 200, bbox_inches = 'tight')
plt.show()







# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% task 5: uniform omegas, fixed thetas, change omegas, t-vs-r repeated %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%time
# --> have same positions as first run, but the natural frequencies get changed

N = 400                  # 2000
T = 200                 # 200
dt = 0.05
runs = 10

numberOfTimes = int(T/dt)
t_range = [round(i * dt, 4) for i in range(numberOfTimes + 1)]
K_range = [1]


fixedThetaPop = OscPopulation('uniform')
# savedThetas = deepcopy(fixedThetaPop)
rLists = []


for run in range(runs):
    _fixedThetaPop = OscPopulation('uniform')

    for n in range(N):
        _fixedThetaPop.list_os[n].currentTheta = fixedThetaPop.list_os[n].currentTheta
        _fixedThetaPop.list_os[n].lastTheta = 0.0
    print('original thetas restored')

    rNewList = _fixedThetaPop.runT()
    rLists.append(rNewList)

    print('############# run', str(run), 'done #############\n')


### graphing ###
plt.figure(figsize = (10,6))
plt.title('Evolution of r')
plt.ylim(0,1)
plt.xlabel('time t')
plt.ylabel('coherence r')

for run in range(runs):
    plt.plot(t_range, rLists[run])

filename = 'graphics/5 t-vs-r' + '_fixedThetas' + '_omegaDistr=' + fixedThetaPop.omegaDistr + '_N=' + str(N) + '.pdf'
plt.savefig(filename, dpi = 200, bbox_inches = 'tight')
plt.show()





# %%
