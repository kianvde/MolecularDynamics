__author__ = 'Kian'
from visual import *
import variables as var
import numpy as np
class VpyAnimate(object):


    def __init__(self, particles, loops):
        self.particles = particles
        self.posMatrix = np.empty((var.numParticles, var.dimension, loops))
        self.speed = 20
        self.sphereRadius = var.boxSize/30.
        self.init_box()
        self.init_balls()





    def init_balls(self):
        self.balls = []
        for i in range(var.numParticles):
            self.balls = self.balls + [sphere(radius=self.sphereRadius, color=color.red)]
            self.balls[i].pos = vector(self.particles.initposs[i,0],self.particles.initposs[i,1],self.particles.initposs[i,2])
            #balls[i].velocity = vector(self.particles.initvelocc[i,0],self.particles.initvelocc[i,1],self.particles.initvelocc[i,2])
    def init_box(self):
        bSize = var.boxSize+self.sphereRadius
        wallR = box (pos=(var.boxSize/2.,var.boxSize/2.,var.boxSize/2.), size=(bSize,bSize,bSize), color = color.blue, opacity = 0.3)


    def buildCoords(self, loopcount, positions):
         self.posMatrix[:, :, loopcount] = positions

    def plot_anim(self):
        self.x = self.posMatrix[:, 0, :] # Npoint , dimension, frame
        self.y = self.posMatrix[:, 1, :]
        self.z = self.posMatrix[:, 2, :]
        while 1:
            for i in range(np.shape(self.x)[1]):
                rate(self.speed)
                for j in range(var.numParticles):
                    self.balls[j].pos = vector(self.x[j,i],self.y[j,i],self.z[j,i])

    def plot_animm(self, positions):
        self.x = positions[:,0]
        self.y = positions[:,1]
        self.z = positions[:,2]
        for j in range(var.numParticles):
            self.balls[j].pos = vector(self.x[j],self.y[j],self.z[j])



