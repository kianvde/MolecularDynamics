__author__ = 'Kian'
from visual import *
import variables as var
import numpy as np
class VpyAnimate(object):

    def __init__(self, particles, loops):
        self.particles = particles
        self.lineData = np.empty((var.numParticles, var.dimension, loops))
        self.init_balls()
    # wallR = box (pos=(0,0,0), size=(var.boxSize,var.boxSize,var.boxSize), color = color.green)




    def init_balls(self):
        balls = []
        for i in range(var.numParticles):
            balls = balls + [sphere(radius=var.boxSize/30., color=color.orange)]
            balls[i].pos = vector(self.particles.initposs[i,0],self.particles.initposs[i,1],self.particles.initposs[i,2])
            balls[i].velocity = vector(self.particles.initvelocc[i,0],self.particles.initvelocc[i,1],self.particles.initvelocc[i,2])

    def buildCoords(self, loopcount, positions):
         self.lineData[:, :, loopcount] = positions
         self.x = self.lineData[:, 0, :]
         self.y = self.lineData[:, 1, :]
         self.z = self.lineData[:, 2, :]

    #def plot_anim(self):
     #   for i in np.shape(self.x)[1]:

