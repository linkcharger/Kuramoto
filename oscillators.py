import math
import random as rand
import numpy as np
import matplotlib.pyplot as plt


T = 50          # 100
h = 0.01
N = 10           # 1000



class Oscillator:

    def __init__(self, initPhase, natFreq):
        self.omega = natFreq
        self.theta = [initPhase]




class OscPopulation:

    def __init__(self):
        # list of OBJECTS!   
        self.list_os = []

        for n in range(N):
            initPhase = rand.uniform(0, 2 * math.pi)
            natFreq = rand.normalvariate(0,1)

            self.list_os.append(Oscillator(initPhase, natFreq))
            print("oscillator ", n, " created")



    ############ calculate r ############
    def calc_r(self):                                       
        
        # complex sum of all phases 
        sum = 0 + 0j
        for n in range(N):
            sum += math.exp(1j * self.list_os[n].theta[-1])                             # [-1] gives last element in list

        sum = 1/N * sum

        r = abs(sum)

        return r



    # step FOR time t -> use t-1 ####################################
    def stepAll(self, t, K):

        for n in range(N):

            # calculate sum for ONE oscillator n
            sum = 0
            for j in range(N):
                sum += math.sin(self.list_os[j].theta[t-1] - self.list_os[n].theta[t-1])            # radians?

            # theta_dot_t for oscillator n
            theta_dot_t = self.list_os[n].omega + K * sum
        
            # new theta for oscillator n
            theta_t = self.list_os[n].theta[t-1] + h * theta_dot_t
            self.list_os[n].theta.append(theta_t)                                      # going down list of objects, pick object, dial into theta list, append






    def run(self, mode):
        rList= []

        ############################################
        if mode == 'K-vs-r':

            for K in range(0, 5, .5):  

                for t in range(0, T, h):
                    
                    self.stepAll(t, K)

                r_crit = self.calc_r()
                rList.append(r_crit)



            # execute at the very end only
            plt.plot(range(0, 5, .5), rList, 'ko')
            plt.title('Bifurcation diagram')
            plt.xlabel('coupling strength K')
            plt.ylabel('coherence r')
            plt.savefig('K-vs-r__K-' + K + '.png', dpi = 200, bbox_inches = 'tight')



        #####################################
        if mode == 't-vs-r':
            for K in [1,2]:  

                for t in range(0, T, h):   
                    # step all oscillators forward INTO timeperiod t
                    self.stepAll(t,K)

                    # calculate the population's coherence at each time step
                    r = self.calc_r()
                    rList.append(r)


                # execute for both runs, individually
                plt.plot(range(0, T, h), rList, 'ko')
                plt.title('Evolution of r, K = ' + K)
                plt.xlabel('time t')
                plt.ylabel('coherence r')
                plt.savefig('t-vs-r__K-' + K + '.png', dpi = 200, bbox_inches = 'tight')







############################################ execute ############################################
population1 = OscPopulation()
population1.run('K-vs-r')


population2 = OscPopulation()
population2.run('t-vs-r')
