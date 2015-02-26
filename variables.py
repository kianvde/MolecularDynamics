### the variables used in the simulation stored in a single file

## imports
import numpy as np
from potentials import lennardJones
import matplotlib
matplotlib.use("qt4agg")
from matplotlib import pyplot as plt

## Constants and variables are given in nm, ps units (i.e distance 1 = 1 nm, time 1 = 1 ps)

## variables ##
numParticlesAxis = 10   # Number of particles per axis, rescale everything else to accommodate the density
deltaT = 0.032          # time step
density = 0.65          # Density
T = 1.036               # Temperature (in Kelvin)

## constants ##
dimension = 3                               # dimensionality of the system
numParticles = numParticlesAxis**3          # number of particles
boxSize = (numParticles/density)**(1./3)    # length of the box side of the box
m = 48                                      # Mass
kB = 1.0                                    # Boltzmann constant
a = ((kB * T) / m)**0.5                     # Maxwell-Boltzmann standard deviation per component sqrt(3kT/m)
eps =  1.0                                  # Lennard-Jones depth of potential well
sigma = 1.0                                 # Lennard-Jones sigma
rMin = sigma * 2.0**(1.0/6.0)               # Lennard-Jones distance at which potential is minimal
rCutoff = 0.45 * boxSize #3.3*sigma                         # interaction range for the particles
tau = deltaT / 0.25                         # Thermostat rise time


## Particles
class Particles(object):

    ## initialization ##
    def __init__(self):

        self.initPositions()
        self.initVelocities()
        self.force = np.zeros(np.shape(self.positions))
        self.energy = 0
        self.potential = 0
        self.pressure = 0
        self.virialFactor = np.array((0,0))
        self.temperature = T
        self.initposs = self.initPositions()
        self.initvelocc = self.initVelocities()

    # initialize the particle positions
    def initPositions(self):
        # Using cubic lattice
        # volumeBox   = boxSize**dimension

        self.positions = np.zeros((numParticles,dimension))
        numAxis = float(numParticlesAxis)
        increment = int(round(numAxis))
        posAxis     = np.arange(0,numAxis)/numAxis * boxSize
        k=0
        for j in range(0,numParticles,increment):
            self.positions[j:j+increment, 0] = posAxis      # For every n particles that are on an axis, set coordinates of those n. Coords are in posAxis.
                                                            # Let's say these are the x coords, then we have x0,x1..xn,x0,x1...xn etc.
            if (j%increment**2)==0:                         # Here add the 'z' coordinates after n**2
                self.positions[j:j+increment**2, 2] = np.array([posAxis[k]]*increment**2)
                k+=1
                i = 0
            if (j%increment)==0:                  #Add the 'y' coordinates, after n repetitions of x0->xn
                self.positions[j:j+increment, 1] = np.array([posAxis[i]]*increment)
                i += 1
        # A start with using FCC
        # Using fcc lattice, we know density of one fcc cube is 14/a**3 units/m3
        # volumeBox   = boxSize**dimension
        # partDenisty = numParticles/volumeBox
        # sideFcc     = (14./partDenisty)**(1./3)
        return self.positions

    # initialize the particle velocities
    def initVelocities(self):

        # initiate velocities components according to MB distribution for the
        # speed, in nm / ps
        # (i.e. Gaussian distribution with mean=0 and std(=a)=sqrt(3kT/m) for the
        # velocity components
        self.velocities = np.random.normal(0., a, (numParticles,dimension))

    ## update functions ##
    def update(self, dT):

        # # velocity-verlet (2nd order)
        # self.updateParticles(dT, 1)
        # self.updateVelocities(dT, 0.5)
        # self.updateForce()
        # self.updateVelocities(dT, 0.5)

        # 4th-order symplectic RK4 integrator (Forest & Ruth, 1989)
        X = (1.0/6.0) * (2.0**(1.0/3.0) + 2.0**(-1.0/3.0) - 1)
        self.updateParticles(dT, X + 0.5)
        self.updateForce()
        self.updateVelocities(dT, 2.0*X + 1.0)
        self.updateParticles(dT, -X)
        self.updateForce()
        self.updateVelocities(dT, -4.0*X - 1.0)
        self.updateParticles(dT, -X)
        self.updateForce()
        self.updateVelocities(dT, 2.0*X + 1.0)
        self.updateParticles(dT, X + 0.5)
        self.updateForce()

        v2 = np.sum(self.velocities**2,axis=None)
        vavg = (3.0 * kB * T / m)**0.5
        self.temperature = (2.0/3.0) * (0.5 * m * v2) / (numParticles * kB)  # Paper Verlet 1967

        # Rescaling according to Berendsen thermostat
        lambdy = (1.0 + (deltaT / tau) * (T / self.temperature - 1.0))**0.5
        self.velocities = self.velocities * lambdy

        v2 = np.sum(self.velocities**2,axis=None)
        self.energy = (0.5 * m * v2 + self.potential) # Energy in Joule
        self.temperature = (2.0/3.0) * (0.5 * m * v2) / (numParticles * kB)  # Paper Verlet 1967

        self.calculatePressure()

    # update the particle positions
    def updateParticles(self, dT, D):

        self.positions += D * self.velocities * dT #(verlet) + 0.5 * (self.force / m) * (dT**2)

        # translate the particles outside of the box
        # +boxSize if positionComponent < 0, -boxSize if positionComponent > 5
        posTranslation = self.positions < 0
        negTranslation = self.positions > boxSize
        self.positions[posTranslation] += boxSize
        self.positions[negTranslation] -= boxSize

    # update the particle velocities and forces
    def updateVelocities(self, dT, C):

        self.velocities += C * (self.force / m) * dT

    # update the forces on the particles and calculate the potential energy
    def updateForce(self):

        self.force, self.potential, vF = lennardJones(self.positions)
        self.virialFactor = np.append(self.virialFactor, vF)

    # calculate the pressure using the virial theorem
    def calculatePressure(self):

        if (len(self.virialFactor) < 10):
            vF = np.mean(self.virialFactor)
        else:
            vF = np.mean(self.virialFactor[-9:])

        self.pressure = density*T*(kB - vF/(3.0*numParticles*T))

    ## UNUSED ##

    # get the particles within a distance imageSize from the box boundaries and translate them to outside
    # the box for force computation to simulate an infinite volume
    #
    # return -> a by d matrix with all the imaged particles
    def getTranslatedImages(self):

        p = self.positions                          # original particle position matrix
        t = np.zeros((np.shape(p)))                 # matrix containing translation values per coordinate
        t[p < rCutoff] = boxSize
        t[p > boxSize - rCutoff] = -boxSize
        tBool = t != 0
        numWalls = np.sum(tBool, axis=1)             # number of walls the particle is close to

        # select cases with one or more translations and add translations and translate
        # all coordinates close to the wall
        translatedImages = p[numWalls >= 1, :] + t[numWalls >= 1, :]

        # select cases with two or more translations and translate for single
        # coordinates close to the wall
        p2 = p[numWalls >= 2, :]
        t2 = t[numWalls >=2, :]

        for column in range(3):
            pCol = p2[:,column]
            tCol = t2[:,column]
            pConcat = p2[tCol != 0, :]
            pConcat[:, column] += tCol[tCol != 0]
            translatedImages = np.concatenate((translatedImages, pConcat), axis=0)

        # select cases with 3 translations and translate for two
        # coordinates close to the wall
        p3 = p[numWalls == 3, :]
        t3 = t [numWalls == 3, :]
        for column in range(3):
            t2Col = t3.copy()
            t2Col[:, column] = 0
            translatedImages = np.concatenate((translatedImages, p3 + t2Col), axis=0)

        return translatedImages

class plotHelper(object):

    def __init__(self):
        self.fig = plt.figure()
        self.fig2 = plt.figure()
        self.axis = self.fig.add_subplot(111)
        self.axis2 = self.fig2.add_subplot(111)
        self.partTemp = []
        pass

    def plotTemp(self,loopnum, temp):
        self.partTemp = self.partTemp + [temp]
        print self.partTemp
        self.axis.plot(np.arange(loopnum+1), self.partTemp)
        plt.show()

    def plotEnergy(self, loopnum):
        self.axis2.plot(loopnum, Particles.energy)
        plt.pause(0.00001)
