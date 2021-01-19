# %% setup & construction %%
import math
import random as rand
import numpy as np
import matplotlib.pyplot as plt



T = 50        # 100
h = 0.01
N = 50        # 1000
numberOfPoints = int(T/h)
t_range = [round(i * h, 4) for i in range(1, numberOfPoints)]
K_range = np.linspace(0, 5, 25)             # from 0 to 5, give me 20 points total, including endpoint = 5


## when I still used indices -> indexes need be integers ##
# int_t_range = [round(t_range[i] / h) for i in range(numberOfPoints-1)]
# t = 0 is setup -> in first loop t = 1, will grab t-1 = 0 in sin()!




class Oscillator:

    def __init__(self, initPhase, natFreq):
        self.omega = natFreq
        self.lastTheta = 0
        self.currentTheta = initPhase




class OscPopulation:

    def __init__(self):
        # list of objects  
        self.list_os = []

        for n in range(N):
            initPhase = rand.uniform(0, math.tau)
            natFreq = rand.normalvariate(0,1)

            self.list_os.append(Oscillator(initPhase, natFreq))
        print("oscillators created")



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
        
        # step through time -> hand over value
        for n in range(N):
            self.list_os[n].lastTheta = self.list_os[n].currentTheta

        # calculate new thetas
        for n in range(N):
            # calculate differential sum for ONE oscillator n for time t
            sum = 0
            for j in range(N):
                sum += math.sin(self.list_os[j].lastTheta - self.list_os[n].lastTheta)                  

            # theta_dot_t for oscillator n 
            theta_dot_t = self.list_os[n].omega + K/N * sum
        
            # new theta for oscillator n -> euler step
            theta_t = self.list_os[n].lastTheta + h * theta_dot_t

            # going down list of objects, pick object, dial into theta list, append
            self.list_os[n].currentTheta = theta_t






    def run(self, mode):
        

        ####################################################################################################################################
        if mode == 'K-vs-r':
            r_critList= []
                            
            for K in K_range:  

                for t in t_range:
                    self.stepAll(K)

                r_crit = self.calc_r()
                r_critList.append(r_crit)
                print('K = ' + str(round(K,3)) + ' | r_crit = ' + str(round(r_crit, 3)))
                


            # execute at the very end only
            plt.figure(figsize = (10,6))
            plt.title('Bifurcation diagram')
            plt.plot(K_range, r_critList, 'ko')
            plt.xlabel('coupling strength K')
            plt.ylabel('coherence r')
            plt.savefig('K-vs-r_K-' + '_N' + str(N) + '_T' + str(T) + '.pdf', dpi = 200, bbox_inches = 'tight')
            plt.show()



        ####################################################################################################################################
        if mode == 't-vs-r':
            rList = [[], [], []]
            
            for K in [1,2]:  

                for t in t_range:
                    # status
                    if t % 10 == 0: print('t =', t, 'done') 


                    # step all oscillators forward INTO timeperiod t
                    self.stepAll(K)

                    # calculate the population's coherence at each time step
                    r = self.calc_r()
                    rList[K].append(r)
                
                print('K =', K, 'done')
                    

            ### graphing ###

            ## individual plots ## 
            # plt.figure(figsize = (10,12))
            # plt.subplot(2, 1, 1)                                        # 2 rows, 1 column; this is the first
            # plt.plot(t_range, rList[1], 'k', label = 'K = 1')
            # plt.legend()
            # plt.subplot(2, 1, 2)                                        # 2 rows, 1 column; this is the second
            # plt.plot(t_range, rList[2], 'c', label = 'K = 2')
            # plt.legend()

            ## in one plot ##
            plt.figure(figsize = (10,6))
            plt.title('Evolution of r')
            plt.plot(t_range, rList[1], 'k', label = 'K = 1')
            plt.plot(t_range, rList[2], 'c', label = 'K = 2')
            plt.legend()
            
            plt.xlabel('time t')
            plt.ylabel('coherence r')
            plt.savefig('t-vs-r_K-' + '_N' + str(N) + '_T' + str(T) +'.pdf', dpi = 200, bbox_inches = 'tight')
            plt.show()






############################################ execute ############################################





# %% task 2 %%
pop2 = OscPopulation()
pop2.run('t-vs-r')




# %% task 1 %%
pop1 = OscPopulation()
pop1.run('K-vs-r')