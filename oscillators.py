# %% setup & construction %%
import math
import random as rand
import numpy as np
import matplotlib.pyplot as plt



############################################### options for optimisation #############################################
#
# jit
# double sum -> even better, average phase, then ...? wait what?
# 
#####################################################################################################################





class Oscillator:

    def __init__(self, initPhase, omega):
        self.omega = omega
        self.lastTheta = 0
        self.currentTheta = initPhase




class OscPopulation:

    def __init__(self, omegaDistr = 'normal'):
        # list of objects  
        self.list_os = []
        self.omegaDistr = omegaDistr
       
       # complete initialisation
        for n in range(N):
            if omegaDistr == 'normal': omega = rand.normalvariate(0,1)
            elif omegaDistr == 'uniform': omega = rand.uniform(-0.5, 0.5)
                    
            initPhase = rand.uniform(0, math.tau)        
            
            self.list_os.append(Oscillator(initPhase, omega))

        print('oscillators initialised with', omegaDistr, 'distibution of omegas' )



    def changeThetas(self):
        for n in range(N):
            initPhase = rand.uniform(0, math.tau)

            # update thetas
            self.list_os[n].lastTheta = 0
            self.list_os[n].currentTheta = initPhase
            
            # dont touch omegas


    def changeOmegas(self):
        for n in range(N):
            omega = rand.uniform(-0.5, 0.5)

            # update omegas
            self.list_os[n].omega = omega
           
            # dont touch thetas

        



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
        sum = 0 + 0j
        for n in range(N):
            sum += math.e**(1j * self.list_os[n].currentTheta)                             
        sum = sum / N 

        # r = math.sqrt(sum.real**2 + sum.imag**2)
        r = abs(sum)
        
        return r



    # step FOR time t -> use t-1 ####################################
    def stepAll(self, K):
        
        for n in range(N):                                                                  # step through time -> hand over value
            self.list_os[n].lastTheta = self.list_os[n].currentTheta

        # calculate new thetas
        for n in range(N):
            
            sum = 0
            for j in range(N):                                                              # calculate differential sum for ONE oscillator n for time t
                sum += math.sin(self.list_os[j].lastTheta - self.list_os[n].lastTheta)                  

            theta_dot_t = self.list_os[n].omega + K/N * sum                                 # theta_dot_t for oscillator n
            
            theta_t = self.list_os[n].lastTheta + dt * theta_dot_t                          # new theta for oscillator n -> euler step

            self.list_os[n].currentTheta = theta_t                                          # going down list of objects, pick object, dial into theta list, append






    def run(self, mode, run = 0):
        
        if mode == 'K-vs-r':
            r_critList= []
                            
            for K in K_range:  

                for t in t_range:
                    self.stepAll(K)

                r_crit = self.calc_r()
                r_critList.append(r_crit)
                print('K = ' + str(round(K,3)) + ' | r_crit = ' + str(round(r_crit, 3)))
                
            return r_critList



        if mode == 't-vs-r':
                rList = [[]]
                
                for K in K_range:  
                    rList.append([])

                    for t in t_range:
                        # status
                        if t % 10 == 0: print('t =', t, 'done') 


                        # step all oscillators forward INTO timeperiod t
                        self.stepAll(K)

                        # calculate the population's coherence at each time step
                        r = self.calc_r()
                        rList[K].append(r)
                        
                return rList
                





# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% task 1: normal omegas, K-vs-r %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%time


N = 20             # 1000  -  100 took about 45min
T = 100             # 100
dt = 0.01
K = 4
dk = 0.1

numberOfTimes = int(T/dt)
numberOfK = int(K/dk)
t_range = [round(i * dt, 4) for i in range(numberOfTimes + 1)]
K_range =  [round(i * dk, 4) for i in range(numberOfK + 1)]            


pop1 = OscPopulation()
r_critList = pop1.run('K-vs-r') 


# graphics 
plt.figure(figsize = (10,6))
plt.title('System state diagram')
plt.xlim(0, 4)
plt.plot(K_range, r_critList, 'ko')
plt.xlabel('coupling strength K')
plt.ylabel('coherence r')
plt.savefig('graphics/K-vs-r_' + 'omegaDistr=' + pop1.omegaDistr + '_N=' + str(N) + '_T=' + str(T) + '.pdf', dpi = 200, bbox_inches = 'tight')
plt.show()








# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% task 2: normal omegas, t-vs-r %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%time


N = 20             # 1000
T = 100             # 100
dt = 0.01

numberOfTimes = int(T/dt)
t_range = [round(i * dt, 4) for i in range(numberOfTimes + 1)]      
K_range = [1, 2]

pop2 = OscPopulation('normal')
rList = pop2.run('t-vs-r')

### graphing ###
plt.figure(figsize = (10,6))
plt.title('Evolution of r')
plt.plot(t_range, rList[1], 'k', label = 'K = 1')
plt.plot(t_range, rList[2], 'c', label = 'K = 2')
plt.legend()

plt.xlabel('time t')
plt.ylabel('coherence r')
plt.savefig('graphics/t-vs-r_' + 'omegaDistr=' + pop2.omegaDistr + '_N=' + str(N) + '_T=' + str(T) + '.pdf', dpi = 200, bbox_inches = 'tight')
plt.show()








# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% task 3: uniform omegas, K-vs-r %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%time


N = 20                   # 2000    - 100, 100, took 6min
T = 100                   # 200
dt = 0.05
K = 1.5
dk = 0.03

numberOfTimes = int(T/dt)
numberOfK = int(K/dk)
t_range = [round(i * dt, 4) for i in range(numberOfTimes + 1)]
K_range =  [round(i * dk, 4) for i in range(numberOfK + 1)]  

pop3 = OscPopulation('uniform')
r_critList = pop3.run('K-vs-r') 


# graphics 
plt.figure(figsize = (10,6))
plt.title('System state diagram')
plt.xlim(0, 4)
plt.plot(K_range, r_critList, 'ko')
plt.xlabel('coupling strength K')
plt.ylabel('coherence r')
plt.savefig('graphics/K-vs-r_' + 'omegaDistr=' + pop3.omegaDistr + '_N=' + str(N) + '_T=' + str(T) + '.pdf', dpi = 200, bbox_inches = 'tight')
plt.show()






# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% task 4: uniform omegas, fixed omegas, change thetas, t-vs-r repeated %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%time
# --> same natural frequencies as first time, but start at different positions every time

N = 50                      # 2000
T = 50                      # 200
dt = 0.05

numberOfTimes = int(T/dt)
t_range = [round(i * dt, 4) for i in range(numberOfTimes + 1)]
K_range = [1]


fixedOmegaPop = OscPopulation('uniform')
resetPop = fixedOmegaPop
rLists = []

for run in range(10):
    fixedOmegaPop = resetPop
    print('## population reset ##')

    fixedOmegaPop.changeThetas()
    print('## thetas changed ##')

    rNewList = fixedOmegaPop.run('t-vs-r')
    rLists.append(rNewList[1])

    print('############# run', str(run), 'done #############\n')



### graphing ###
plt.figure(figsize = (10,6))
plt.title('Evolution of r')
plt.ylim(0,1)
plt.xlabel('time t')
plt.ylabel('coherence r')

for run in range(10):
    plt.plot(t_range, rLists[run])

filename = 'graphics/t-vs-r' + '_fixedOmegas' + '_omegaDistr=' + fixedOmegaPop.omegaDistr + '_N=' + str(N) + '_T=' + str(T) + '.pdf'
plt.savefig(filename, dpi = 200, bbox_inches = 'tight')
plt.show()









# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% task 5: uniform omegas, fixed thetas, change omegas, t-vs-r repeated %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%time
# --> have same positions as first run, but the natural frequencies get changed

N = 50                  # 2000
T = 50                 # 200
dt = 0.05

numberOfTimes = int(T/dt)
t_range = [round(i * dt, 4) for i in range(numberOfTimes + 1)]
K_range = [1]


fixedThetaPop = OscPopulation('uniform')
resetPop = fixedThetaPop
rLists = []


for run in range(10):
    fixedThetaPop = resetPop
    print('## population reset ##')

    fixedThetaPop.changeOmegas()
    print('## omegas changed ##')

    rNewList = fixedThetaPop.run('t-vs-r')
    rLists.append(rNewList[1])

    print('############# run', str(run), 'done #############\n')


### graphing ###
plt.figure(figsize = (10,6))
plt.title('Evolution of r')
plt.ylim(0,1)
plt.xlabel('time t')
plt.ylabel('coherence r')

for run in range(10):
    plt.plot(t_range, rLists[run])

filename = 'graphics/t-vs-r' + '_fixedThetas' + '_omegaDistr=' + fixedThetaPop.omegaDistr + '_N=' + str(N) + '_T=' + str(T) + '.pdf'
plt.savefig(filename, dpi = 200, bbox_inches = 'tight')
plt.show()


            