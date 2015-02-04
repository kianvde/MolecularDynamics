import numpy as np
import variables as var

class MolDynSimulation(object) :

    # move to parameter file
    numParticles = 1000
    dimension = 2

    # position and velocity vectors for the particles
    particlePos = np.zeros([var.numParticles, var.dimension])
    particleVel = np.zeros([var.numParticles, var.dimension])

    # all the initialization functions are called in this block
    def __init__(self) :
        pass


    # start of the simulation
    def start(self):
        while (True):
            pass



# init and loop
molDynSimulation = MolDynSimulation()
molDynSimulation.start()
