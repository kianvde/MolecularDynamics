import variables as var
import numpy as np
import animate as an

class MolDynSimulation(object) :

    # all the initialization functions are called in this block
    def __init__(self) :
        self.particles = var.Particles()
        pass


    # start of the simulation
    def start(self):
        np.set_printoptions(precision = 2)
        print(self.particles.positions)
        # Number of simulation loops
        loops = 2
        #Create animation object
        self.animation = an.Animate(var.boxSize, loops, var.dimension, var.numParticles)

        for i in range(loops):
            #Build coordinate matrix for every iteration of the loop
            self.animation.buildLines(i, self.particles.positions)
            # update particles
            self.particles.update(var.deltaT)
        #Resize axis and do animation
        self.animation.plot_anim()


# init and loop
molDynSimulation = MolDynSimulation()
molDynSimulation.start()
