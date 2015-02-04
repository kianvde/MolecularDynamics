import Particles

class MolDynSimulation(object) :

    particles = Particles
    # all the initialization functions are called in this block
    def __init__(self) :
        self.particles = Particles()
        pass


    # start of the simulation
    def start(self):

        # simulation loop
        while (True):

            # update particles
            self.particles.update()
            pass



# init and loop
molDynSimulation = MolDynSimulation()
molDynSimulation.start()
