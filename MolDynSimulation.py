import variables as var

class MolDynSimulation(object) :

    # all the initialization functions are called in this block
    def __init__(self) :
        self.particles = var.Particles()
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
