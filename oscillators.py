# %% setup & construction %%
import math
import random as rand
import numpy as np
import matplotlib.pyplot as plt
from numba import jit, double
import time
from scipy.integrate import quad


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

        print('## oscillators initialised with', omegaDistr.upper(), 'distibution of omegas ##')




    ############ calculate r ############
    def calc_r(self):                                       
        
        ###### one way to do it #######
        # complex sum of all phases 
        realSum = 0.0 
        imaginarySum = 0.0

        for n in range(N): 
            realSum += math.cos(self.list_os[n].currentTheta)/N
            imaginarySum += math.sin(self.list_os[n].currentTheta)/N
        
        # realSum = realSum/N 
        # imaginarySum = imaginarySum/N
        r = math.sqrt(realSum * realSum + imaginarySum * imaginarySum)


        ###### another way to do it #######
        # sum = 0.0 + 0.0j

        # for n in range(N): sum += math.e**(1.0j * self.list_os[n].currentTheta)                             
        # sum = sum / N 

        # # r = math.sqrt(sum.real**2 + sum.imag**2)
        # r = abs(sum)
        
        return r

    

    # wrapper method
    def oneStepForAll(self, K):

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

        # calculating
        newLastThetaArray, newCurrentThetaArray = self._oneStepForAll(lastThetaArray, currentThetaArray, omegaArray, K)

        # setter: back to OOP
        for n in range(N):
            self.list_os[n].lastTheta = newLastThetaArray[n]
            self.list_os[n].currentTheta = newCurrentThetaArray[n]

    
    
    # step FOR time t -> use t-1 ####################################
    @staticmethod
    @jit(nopython = True)            #('double[:](double[:], double[:], double[:], int8)', nopython = True)
    def _oneStepForAll(lastThetaArray, currentThetaArray, omegaArray, K):
        
        for n in range(N):                                                                  # step through time -> hand over value
            lastThetaArray[n] = currentThetaArray[n]
            currentThetaArray[n] = 0.0

        # calculate new thetas
        for n in range(N):
            
            sum = 0.0
            for j in range(N): sum += math.sin(lastThetaArray[j] - lastThetaArray[n])                  

            theta_dot_t = omegaArray[n] + K/N * sum                                 # theta_dot_t for oscillator n
            
            theta_t = lastThetaArray[n] + dt * theta_dot_t                          # new theta for oscillator n -> euler step

            currentThetaArray[n] = theta_t                                          # going down list of objects, pick object, dial into theta list, append

        return lastThetaArray, currentThetaArray




    def runK(self):
        
        _r_critList= []
                        
        for K in K_range:  

            # re-randomise
            for n in range(N):
                # finding new random values
                # omega
                if self.omegaDistr == 'normal': newOmega = rand.normalvariate(0.0,1.0)
                elif self.omegaDistr == 'uniform': newOmega = rand.uniform(-0.5, 0.5)  
                # theta
                newInitPhase = rand.uniform(0, math.tau)      
                
                # setting them
                self.list_os[n].omega = newOmega
                self.list_os[n].lastTheta = 0.0                                                 # gets value immediately, but instantiate type
                self.list_os[n].currentTheta = newInitPhase                                    # handed over immediately on first run to lastTheta
            # print('## 3 variables re-randomised ##' )  

            # run
            for t in t_range: self.oneStepForAll(K)

            # calc r
            _r_crit = self.calc_r()
            _r_critList.append(_r_crit)
            print('K = ' + str(round(K,3)) + ' | _r_crit = ' + str(round(_r_crit, 3)))
            
        return _r_critList
                





    def runT(self, K):
        _rList = []

        for t in t_range:
            self.oneStepForAll(K)

            r = self.calc_r()                               # calculate the population's coherence at each time step
            _rList.append(r)
                    
            # status
            if t % 10 == 0: print('t =', t, 'done') 
    
    
        return _rList
                




def g(w):
    return 1/(math.sqrt(2 * math.pi)) * math.e**(-1/2 * w**2)




def calc_r_integ(K):
    r = 1 - dr 
    RHS = 0

    while r > 0:
        
        def integrand(theta): 
            return (math.cos(theta))**2 * g(K * r * math.sin(theta))

        RHS = K * quad(integrand, -(math.pi/2), math.pi/2)[0] 
        # print(r, RHS)

        if RHS > 0 and (RHS < 1 and RHS > (1 - rHysteresis)) or (RHS > 1 and RHS < (1 - rHysteresis)):
            break
        else: 
            r -= dr

    if K % 0.5 == 0: print('K =', K, '| final RHS =', RHS, ' | r =', r)


    return r









# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% task 3: uniform omegas, K-vs-r %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%time
N = 2000                   # 2000    - 100, 100, took 6min
T = 200                   # 200
dt = 0.05
Kmax = 1.5
dk = 0.03

numberOfTimes = int(T/dt)
numberOfK = int(Kmax/dk)
t_range = [round(i * dt, 4) for i in range(numberOfTimes + 1)]
K_range =   [round(i * dk, 4) for i in range(numberOfK + 1)]  # np.linspace(0, Kmax, numberOfK + 1) #
# K_range = [0.5, 0.7, 1.0, 1.2] 

pop3 = OscPopulation('uniform')
r_critList = pop3.runK() 


#################################### graphing #################################### 
plt.figure(figsize = (10,6))
plt.title('System state diagram')
plt.xlim(0, 4)
plt.ylim(0, 1)
plt.xlabel('coupling strength K')
plt.ylabel('coherence r')

plt.plot(K_range, r_critList, 'ko')
filename = 'graphics/3_K-vs-r' + '_omegaDistr=' + pop3.omegaDistr + '_N=' + str(N) + '_' + str(int(time.time())) + '.pdf'
plt.savefig(filename, dpi = 200, bbox_inches = 'tight')
plt.show()


#################################### graphing #################################### 
plt.figure(figsize = (10,6))
plt.title('System state diagram')
plt.ylim(0, 1)
plt.xlabel('coupling strength K')
plt.ylabel('coherence r')

plt.plot(K_range, r_critList, 'ko')
filename = 'graphics/3_K-vs-r_' + 'omegaDistr=' + pop3.omegaDistr + '_zoomed' + '_N=' + str(N) + '_' + str(int(time.time())) + '.pdf'
plt.savefig(filename, dpi = 200, bbox_inches = 'tight')
plt.show()





# %%%%%%%%%%%%%%%%% task 4: uniform omegas, fixed omegas, change thetas, t-vs-r repeated %%%%%%%%%%%%%%%%%
#%time
# --> same natural frequencies as first time, but start at different positions every time

N = 2000                      # 2000
T = 200                      # 200
dt = 0.05
runs = 10

numberOfTimes = int(T/dt)
t_range = [round(i * dt, 4) for i in range(numberOfTimes + 1)] # np.linspace(0, T, numberOfTimes + 1)
K = 1             # 1

fixedOmegaPop = OscPopulation('uniform')
rLists = []


for run in range(runs):
    tempPop = OscPopulation('uniform')

    for n in range(N): tempPop.list_os[n].omega = fixedOmegaPop.list_os[n].omega
    print('original omegas restored')

    rNewList = tempPop.runT(K)
    rLists.append(rNewList)

    print('############# run', str(run), 'done #############\n') 


#################################### graphing ####################################
plt.figure(figsize = (10,6))
plt.title('Evolution of r')
plt.ylim(0,1)
plt.xlabel('time t')
plt.ylabel('coherence r')

for run in range(runs):
    # for t in t_range: rLists[run][int(t/dt)] = 2 * rLists[run][int(t/dt)]
    plt.plot(t_range, rLists[run])

filename = 'graphics/4_t-vs-r' + '_fixedOmegas' + '_omegaDistr=' + fixedOmegaPop.omegaDistr + '_N=' + str(N) + '_' + str(int(time.time())) + '.pdf'
plt.savefig(filename, dpi = 200, bbox_inches = 'tight')
plt.show()







# %%%%%%%%%%%%%%%%% task 5: uniform omegas, fixed thetas, change omegas, t-vs-r repeated %%%%%%%%%%%%%%%%%
#%time
# --> have same positions as first run, but the natural frequencies get changed

N = 2000                  # 2000
T = 200                 # 200
dt = 0.05
runs = 10

numberOfTimes = int(T/dt)
t_range = np.linspace(0, T, numberOfTimes + 1) #[round(i * dt, 4) for i in range(numberOfTimes + 1)]
K = 1


fixedThetaPop = OscPopulation('uniform')
# savedThetas = deepcopy(fixedThetaPop)
rLists = []


for run in range(runs):
    _fixedThetaPop = OscPopulation('uniform')

    for n in range(N):
        _fixedThetaPop.list_os[n].currentTheta = fixedThetaPop.list_os[n].currentTheta
        _fixedThetaPop.list_os[n].lastTheta = 0.0
    print('original thetas restored')

    rNewList = _fixedThetaPop.runT(K)
    rLists.append(rNewList)

    print('############# run', str(run), 'done #############\n')


#################################### graphing ####################################
plt.figure(figsize = (10,6))
plt.title('Evolution of r')
plt.ylim(0,1)
plt.xlabel('time t')
plt.ylabel('coherence r')

for run in range(runs):
    plt.plot(t_range, rLists[run])

filename = 'graphics/5_t-vs-r' + '_fixedThetas' + '_omegaDistr=' + fixedThetaPop.omegaDistr + '_N=' + str(N) + '_' + str(int(time.time())) + '.pdf'
plt.savefig(filename, dpi = 200, bbox_inches = 'tight')
plt.show()








# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% task 1: normal omegas, K-vs-r %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%time

#################################### simulation ####################################
N = 1000             # 1000 
T = 100             # 100
dt = 0.01
Kmax = 4
dk = 0.1

numberOfTimes = int(T/dt)
numberOfK = int(Kmax/dk)
t_range = [round(i * dt, 4) for i in range(numberOfTimes + 1)]
K_range =  [round(i * dk, 4) for i in range(numberOfK + 1)]            


pop1 = OscPopulation('normal')
rListSim = pop1.runK() 


#################################### integration ####################################
dk = 0.01
numberOfK = int(Kmax/dk)
K_rangeInteg =  [round(i * dk, 4) for i in range(numberOfK + 1)]  

rHysteresis = 0.01
dr = 0.001


rListInteg = []
for K in K_rangeInteg:
    r_crit = calc_r_integ(K)
    rListInteg.append(r_crit)
    # print('K = ' + str(K) + ' | r = ' + str(r_crit))


#################################### graphing ####################################
plt.figure(figsize = (10,6))
plt.title('System state diagram')
plt.xlim(0, 4)
plt.ylim(0, 1)
plt.xlabel('coupling strength K')
plt.ylabel('coherence r')

plt.plot(K_range, rListSim, 'ko', label = 'simulated')
plt.plot(K_rangeInteg, rListInteg, label = 'integrated')
plt.legend()
filename = 'graphics/1_K-vs-r' + '_omegaDistr=' + pop1.omegaDistr + '_N=' + str(N) + '_' + str(int(time.time())) + '.pdf'
plt.savefig(filename, dpi = 200, bbox_inches = 'tight')
plt.show()







# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% task 2: normal omegas, t-vs-r %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%time

N = 1000             # 1000
T = 100             # 100
dt = 0.01
K_range = [1, 2, 3]

numberOfTimes = int(T/dt)
t_range = [round(i * dt, 4) for i in range(numberOfTimes + 1)]      

rLists = [[]]

for K in K_range:
    pop2 = OscPopulation('normal')
    rList = pop2.runT(K)

    rLists.append(rList)





#################################### graphing ####################################
plt.figure(figsize = (10,6))
plt.title('Evolution of r')
plt.ylim(0, 1)
plt.xlabel('time t')
plt.ylabel('coherence r')

for K in K_range:
    plt.plot(t_range, rLists[K], label = 'K = ' + str(K))
plt.legend()
filename = 'graphics/2_t-vs-r' + '_omegaDistr=' + pop2.omegaDistr + '_N=' + str(N) + '_' + str(int(time.time())) + '.pdf'
plt.savefig(filename, dpi = 200, bbox_inches = 'tight')
plt.show()








