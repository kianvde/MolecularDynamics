import variables as var
import numpy as np
import vpy_animate as an

class MolDynSimulation(object) :

    # all the initialization functions are called in this block
    def __init__(self) :
        self.particles = var.Particles()

        # Number of simulation loops
        self.numIterations = 1000

        # Create animation object
        #self.animation = an.Animate(var.boxSize, self.numIterations, var.dimension, var.numParticles)
        self.animation = an.VpyAnimate(self.particles, self.numIterations)
        self.plotGraph = var.plotHelper()


    # start of the simulation
    def start(self):

        for i in range(self.numIterations):
            # Build coordinate matrix for every iteration of the loop
            #self.animation.buildCoords(i, self.particles.positions)
            # update particles
            self.plotGraph.plotTemp(i)
            self.particles.update(var.deltaT)
            # Set all particles to new position for every update
            #self.animation.plot_anim(self.particles.positions)
            #print i

        # Resize axis and do animation


# init and loop
molDynSimulation = MolDynSimulation()
molDynSimulation.start()
