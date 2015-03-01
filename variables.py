### the variables used in the simulation stored in a single file

## imports
import numpy as np
from potentials import lennardJones
import matplotlib
matplotlib.use("qt4agg")
from matplotlib import pyplot as plt

## Constants and variables are given in natural units (i.e sigma = 1, kB  = 1, epsilon = 1)

## variables ##
numParticlesAxis = 10   # Number of particles per axis
dt = 0.032              # time step
density = 0.65          # Density
T = 1.0                 # Temperature

## constants ##
dimension = 3                               # dimensionality of the system
N = numParticlesAxis**3                     # number of particles
boxSize = (N/density)**(1./3)               # length of the box side of the box
m = 48                                      # Mass
k = 1.0                                     # Boltzmann constant
a = ((k * T) / m)**0.5                      # Maxwell-Boltzmann standard deviation per component sqrt(kT/m)
eps =  1.0                                  # Lennard-Jones depth of potential well
sig = 1.0                                   # Lennard-Jones sigma
rMin = sig * 2.0**(1.0/6.0)                 # Lennard-Jones distance at which potential is minimal
rC = 0.45 * boxSize                         # interaction range for the particles
tau = dt / 0.25                             # Thermostat rise time
numBins = 100                               # number of bins for correlation function calculation

def calculateDerivedConstants():
    global N, boxSize, a, rC, tau
    N = numParticlesAxis**3                     # number of particles
    boxSize = (N/density)**(1./3)               # length of the box side of the box
    a = ((k * T) / m)**0.5                      # Maxwell-Boltzmann standard deviation per component sqrt(kT/m)
    rC = 0.45 * boxSize                         # interaction range for the particles
    tau = dt / 0.25                             # Thermostat rise time


## Particles
# class with the simulation information and update functionality
class Particles(object):

    ## initialization ##
    def __init__(self):

        self.positions = np.zeros((N,dimension))    # storing the particle positions
        self.velocities = np.zeros((N,dimension))   # storing the particle velocity components

        self.force = np.zeros(np.shape(self.positions))
        self.energy = 0
        self.kinetic = 0
        self.potential = 0
        self.pressure = 0
        self.temperature = T

        self.bin = np.zeros(numBins)            # bin for correlation function calculation
        self.virialFactor = np.array((0,0))     # vector for storing virial theorem factors

        self.initPositions()
        self.initVelocities()

    # initialize the particle positions using a cubic lattice
    def initPositions(self):

        numAxis = float(numParticlesAxis)
        increment = int(round(numAxis))
        posAxis     = np.arange(0,numAxis)/numAxis * boxSize
        k=0
        for j in range(0,N,increment):
            self.positions[j:j+increment, 0] = posAxis      # For every n particles that are on an axis, set coordinates of those n. Coords are in posAxis.
                                                            # Let's say these are the x coords, then we have x0,x1..xn,x0,x1...xn etc.
            if (j%increment**2)==0:                         # Here add the 'z' coordinates after n**2
                self.positions[j:j+increment**2, 2] = np.array([posAxis[k]]*increment**2)
                k+=1
                i = 0
            if (j%increment)==0:                  #Add the 'y' coordinates, after n repetitions of x0->xn
                self.positions[j:j+increment, 1] = np.array([posAxis[i]]*increment)
                i += 1

    # initiate velocities components according to MB distribution
    # for the speed
    def initVelocities(self):
        self.velocities = np.random.normal(0., a, (N,dimension))
        self.initialVelocities = self.velocities.copy()


    ## update functions ##
    def update(self, dT, rescale):

        # # velocity-verlet (2nd order)
        # self.updateParticles(dT, 1)
        # self.updateVelocities(dT, 0.5)
        # self.updateForce(True)
        # self.updateVelocities(dT, 0.5)

        # symplectic RK4 integrator (4th order) (Forest & Ruth, 1989)
        X = (1.0/6.0) * (2.0**(1.0/3.0) + 2.0**(-1.0/3.0) - 1)
        self.updateParticles(dT, X + 0.5)
        self.updateForce(False)
        self.updateVelocities(dT, 2.0*X + 1.0)
        self.updateParticles(dT, -X)
        self.updateForce(False)
        self.updateVelocities(dT, -4.0*X - 1.0)
        self.updateParticles(dT, -X)
        self.updateForce(False)
        self.updateVelocities(dT, 2.0*X + 1.0)
        self.updateParticles(dT, X + 0.5)
        self.updateForce(True)

        # Temperature rescaling according to Berendsen thermostat
        if rescale == 1:
            lambdy = (1.0 + (dT/tau)*(T/self.temperature - 1.0))**0.5
            self.velocities = self.velocities * lambdy

        # calculation energy and temperature
        v2 = np.sum(self.velocities**2)
        self.kinetic = 0.5*m*v2
        self.energy = self.kinetic + self.potential
        self.temperature = (2.0/3.0)*(0.5*m*v2)/(N*k)

        # calculate pressure
        self.calculatePressure()

    # update the particle positions
    def updateParticles(self, dT, vFactor):

        self.positions += vFactor*self.velocities*dT

        # for verlet
        if (vFactor == 1):
            self.positions += 0.5 * (self.force / m) * (dT**2)

        # translate the particles outside of the box
        self.positions[self.positions < 0] += boxSize
        self.positions[self.positions > boxSize] -= boxSize

    # update the particle velocities and forces
    def updateVelocities(self, dT, C):

        self.velocities += C * (self.force / m) * dT

    # update the forces on the particles and calculate the potential energy
    def updateForce(self, addToCorrelationBin):

        self.force, self.potential, vF, cBin = lennardJones(self.positions)
        # self.virialFactor = np.append(self.virialFactor, vF)
        tailcorrU = 4.0 * np.pi * rC**3 * ((1.0/9.0) * (2.0**(1.0/6.0)/rC)**12 - (1.0/3.0)*(2.0**(1.0/6.0)/rC)**6);
        self.potential = self.potential + tailcorrU

        if (addToCorrelationBin):
            self.bin += cBin
            self.virialFactor = np.append(self.virialFactor, vF)


    ## calculation functions
    # calculate the pressure using the virial theorem
    def calculatePressure(self):

        if (len(self.virialFactor) < 10):
            vF = np.mean(self.virialFactor)
        else:
            vF = np.mean(self.virialFactor[-9:])

        tailcorrfact = - 12.0 * eps * rC
        tailcorr = tailcorrfact * ((1.0 / 11.0) * (rMin / rC)**12 - (1.0 / 5.0) * (rMin / rC)**6)
        self.pressure = density*T*(k - vF/(3.0*N*T)) - density * tailcorr / 3.0

    # get the velocity correlation
    def getVelCorrelation(self):
        return np.sum(self.initialVelocities*self.velocities, axis=1)

    # return the correlation function using the binned distances collected from fortran
    def getCorrelationFunction(self, numIterations):
        i = float(2*numBins)
        r2 = np.linspace(rC/i,rC*((i-1.0)/i),num=numBins)**2

        prefactor = 2.0/(density*(N-1.0)) * 1.0/(4*np.pi*(rC/numBins))

        return (prefactor/(numIterations)) * (self.bin/r2)
