import variables as var
import numpy as np
import vpy_animate as an
import matplotlib.pyplot as plt

class MolDynSimulation(object) :

    # all the initialization functions are called in this block
    def __init__(self) :
        self.particles = var.Particles()

        # Number of simulation loops
        self.numIterations = 1000
        self.Eold = float("Inf")
        self.Ediff = float("Inf")

        # Create animation object
        # self.animation = an.VpyAnimate(self.particles, 1000)


    # start of the simulation
    def start(self):

        energy = np.empty(self.numIterations)
        cV = np.empty(self.numIterations)
        temperature = np.empty(self.numIterations)
        compressibility = np.empty(self.numIterations)
        potentialEnergy = np.empty(self.numIterations)

        e2Avg = 0.0
        eAvg = 0.0
        for i in range(0,self.numIterations):
            self.particles.update(var.deltaT)

            energy[i] = self.particles.energy
            temperature[i] = self.particles.temperature
            potentialEnergy[i] = self.particles.potential
            compressibility[i] = self.particles.pressure/(var.T*var.density)

            cvIndex = i - self.numIterations/2
            if (i >= self.numIterations/2):
                cV[cvIndex] = np.var(energy[self.numIterations/2-10:i])/(var.T**2)

            # if (i<1000):
            #      self.animation.plot_anim(self.particles.positions)

            print i

        plt.plot(compressibility)
        plt.title("compressibility")
        plt.show()
        plt.plot(temperature)
        plt.title("temperature")
        plt.show()
        plt.plot(cV)
        plt.title("cV")
        plt.show()
        plt.plot(energy)
        plt.title("energy")
        plt.show()
        plt.plot(potentialEnergy)
        plt.title("potential energy")
        plt.show()

# init and loop
molDynSimulation = MolDynSimulation()
molDynSimulation.start()
