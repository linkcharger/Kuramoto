import math
import random as rand


T = 100
h = 0.01
N = 1000



class Oscillator:

    def __init__(self, initPhase, natFreq):
        self.natFreq = natFreq
        self.phase = [initPhase]

    def step(self, t):
        sum = 0
        for j in range(N):
            sum += math.sin(oss[j].phase[t] - self.phase[t]) 
            # radians! 
            # time of theta may be t or t+1

        return self.natFreq + sum



def createOscillators():
        
    osses = []
    for n in range(N):
        initPhase = rand.uniform(0, 2 * math.pi)
        natFreq = rand.normalvariate(0,1)

        osses.append(Oscillator(initPhase, natFreq))

    return osses






oss = createOscillators()

def runSim():

    for K in range(0, 5, .5):   

        globalTimer = []
        for t in range (0, T, h):
            globalTimer.append(t)    


            for i in range(N):

                newTheta = oss[i].phase[t-1] + h * oss[i].step(t-1)
                oss[i].phase.append(newTheta)
