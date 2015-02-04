import numpy

class MolDynSimulation(object) :

    # move to parameter file
    numParticles = 1000
    dimension = 2

    # position and velocity vectors for the particles
    particlePos = numpy.zeros([numParticles,dimension])
    particleVel = numpy.zeros([numParticles,dimension])

    # all the initialization functions are called in this block
    def __init__(self) :
        pass


    # start of the simulation
    def start(self):
        while (true):
            pass



# init and loop
molDynSimulation = MolDynSimulation()
molDynSimulation.start()