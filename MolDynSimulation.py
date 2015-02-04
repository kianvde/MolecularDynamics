import variables as var
import numpy as np

class MolDynSimulation(object) :

    # all the initialization functions are called in this block
    def __init__(self) :
        self.particles = var.Particles()
        pass


    # start of the simulation
    def start(self):

        # simulation loop
        for i in range(0, 15):

            # update particles
            self.particles.update(var.deltaT)


# init and loop
molDynSimulation = MolDynSimulation()
molDynSimulation.start()
