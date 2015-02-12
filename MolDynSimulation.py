import variables as var
import numpy as np
import vpy_animate as an

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
        loops = 1000
        #Create animation object
        #self.animation = an.Animate(var.boxSize, loops, var.dimension, var.numParticles)
        self.animation = an.VpyAnimate(self.particles, loops)

        for i in range(loops):
            #Build coordinate matrix for every iteration of the loop
            self.animation.buildCoords(i, self.particles.positions)
            # update particles
            self.particles.update(var.deltaT)
        #Resize axis and do animation
        self.animation.plot_anim()
        #self.particles.

# init and loop
molDynSimulation = MolDynSimulation()
molDynSimulation.start()
