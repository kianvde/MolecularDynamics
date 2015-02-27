import variables as var
import numpy as np
import vpy_animate as an
import matplotlib.pyplot as plt
import time
import os

class MolDynSimulation(object) :

    # all the initialization is done in this block
    def __init__(self) :
        self.particles = var.Particles()

        self.n = 500                                # number of iterations
        self.energy = np.zeros(self.n)              # energy values
        self.cv = np.zeros(self.n)                  # specific heat values
        self.T = np.zeros(self.n)                   # temperature values
        self.potential = np.zeros(self.n)           # potential energy values
        self.compressibility = np.zeros(self.n)     # compressibility values
        self.vCorrelationMean = np.zeros(self.n)
        self.D = np.zeros(self.n)

        # self.animation = an.VpyAnimate(self.particles, self.n)

    # start the simulation loop
    def start(self):

        for i in range(0,self.n):
            self.particles.update(var.dt)

            self.energy[i] = self.particles.energy
            self.potential[i] = self.particles.potential
            self.T[i] = self.particles.temperature
            self.compressibility[i] = self.particles.pressure/(var.T*var.density)

            cvIndex = i - self.n/2
            if (i >= self.n/2):
                self.cv[cvIndex] = np.var(self.energy[(self.n/2-10):i])/(var.T**2)

            self.vCorrelationMean[i] = np.mean(self.particles.getVelCorrelation())
            self.D[i] = np.sum(self.vCorrelationMean * var.dt)

            #      self.animation.plot_anim(self.particles.positions)

            print i

    def write(self):

        if not os.path.exists("data"):
            os.makedirs("data")

        timeString = time.strftime("%m_%d_%H_%M")

        logFile = open('data/log.txt', 'a')
        logFile.write("\n------%s.txt------\n" % timeString)
        logFile.write("##INPUT##\n")
        logFile.write("number of particles (N): %i\n" % var.N)
        logFile.write("density: %g\n" % var.density)
        logFile.write("time step (dt): %g\n" % var.dt)
        logFile.write("initial temperature (T): %g\n" % var.T)
        logFile.write("number of iterations (n): %i\n" % self.n)
        logFile.write("##OUTPUT##\n")
        logFile.write("data file structure from first to last line:\n")
        logFile.write("temperature, energy, potential, compressibility, ")
        logFile.write("specific heat, velocity correlation mean, ")
        logFile.write("diffusion constant, correlation function\n")

        # write the data to a file
        dataFile = open('data/%s.txt' % timeString,'w')
        for value in self.T:
            dataFile.write("%e " % value)
        dataFile.write("\n")
        for value in self.energy:
            dataFile.write("%e " % value)
        dataFile.write("\n")
        for value in self.potential:
            dataFile.write("%e " % value)
        dataFile.write("\n")
        for value in self.compressibility:
            dataFile.write("%e " % value)
        dataFile.write("\n")
        for value in self.cv:
            dataFile.write("%e " % value)
        dataFile.write("\n")
        for value in self.vCorrelationMean:
            dataFile.write("%e " % value)
        dataFile.write("\n")
        for value in self.D:
            dataFile.write("%e " % value)
        dataFile.write("\n")
        for value in self.particles.getCorrelationFunction(self.n):
            dataFile.write("%e " % value)
        dataFile.write("\n")
        dataFile.close()

    # plot the results using matPlotLib
    def plot(self):
        plt.figure()
        plt.plot(self.particles.getCorrelationFunction(self.n))
        plt.title("correlation function")
        plt.figure()
        plt.plot(self.beta)
        plt.title("compressibility")
        plt.figure()
        plt.plot(self.T)
        plt.title("temperature")
        plt.figure()
        plt.plot(self.cv)
        plt.title("cV")
        plt.figure()
        plt.plot(self.energy)
        plt.title("energy")
        plt.figure()
        plt.plot(self.potential)
        plt.title("potential energy")
        plt.figure()
        plt.plot(self.vCorrelationMean)
        plt.title("velocity autocorrelation")
        plt.figure()
        plt.plot(self.D)
        plt.title("diffusion constant")
        plt.show()


# init, loop and plot
molDynSimulation = MolDynSimulation()
molDynSimulation.start()
molDynSimulation.write()
