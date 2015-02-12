### the variables used in the simulation stored in a single file

## imports
import numpy as np
import Potentials as Pot

## constants
## Constants are given in nm, ps units (i.e distance 1 = 1 nm, time 1 = 1 ps)
# dimensionality of the system
dimension = 3

# number of particles in the system
numParticlesAxis = 3
numParticles = numParticlesAxis**3

# time step
deltaT = .1 # 1 microsecond, time step in ps, rescale the rest to fit

# length of the box side of the box
boxSize = 100.0 # Box size in nm, rescale everything else to fit

# part of the box used for particle imaging improving field
imageSize = 0.2*boxSize

# Temperature (in Kelvin)
T = 1.

# Mass
m = 6.64648*10**(-27) # 6.64648*10**(-27) kg

# Boltzmann constant
kB = 1.3806488*10**(-29) # kB = 1.38*10^-29 [nm^2 kg ps^-2 K^-1]

# Maxwell-Boltzmann standard deviation per component sqrt(3kT/m)
a = ((kB * T) / m)**0.5

# Lennard-Jones depth of potential well
eps =  10.22 * kB # Helium Cyrogenics - Steven van Sciver, eps/kB = 10.22 ADDED FACTOR FOR INCREASED INTERACTION
# eps = [J] = [kg m^2 s^-2] = 10^-6 [kg nm^2 ps^-2]

# Lennard-Jones distance at which potential is minimal
rMin = 0.2869 # Helium Cyrogenics - Steven van Sciver, rMin = 0.2869 nm

## classes
class Particles(object):


    ## initialization ##
    def __init__(self):

        self.initposs = self.initPositions()
        self.initvelocc = self.initVelocities()
        self.force = np.zeros(np.shape(self.positions))
        self.energy = 0
        self.potential = 0
        self.temperature = T

    # initialize the particle positions
    def initPositions(self):
        #Using cubic lattice
        #volumeBox   = boxSize**dimension

        self.positions = np.zeros((numParticles,dimension))
        numAxis = float(numParticlesAxis)
        increment = int(round(numAxis))
        posAxis     = np.arange(0,numAxis)/numAxis * boxSize
        k=0
        for j in range(0,numParticles,increment):
            self.positions[j:j+increment, 0] = posAxis      #For every n particles that are on an axis, set coordinates of those n. Coords are in posAxis.
                                                            #Let's say these are the x coords, then we have x0,x1..xn,x0,x1...xn etc.
            if (j%increment**2)==0:                         #Here add the 'z' coordinates after n**2
                self.positions[j:j+increment**2, 2] = np.array([posAxis[k]]*increment**2)
                k+=1
                i = 0
            if (j%increment)==0:                  #Add the 'y' coordinates, after n repetitions of x0->xn
                self.positions[j:j+increment, 1] = np.array([posAxis[i]]*increment)
                i += 1
        #Using fcc lattice, we know density of one fcc cube is 14/a**3 units/m3
        #volumeBox   = boxSize**dimension
        #partDenisty = numParticles/volumeBox
        #sideFcc     = (14./partDenisty)**(1./3)
        return self.positions

    # initialize the particle velocities
    def initVelocities(self):

        # initiate velocities components according to MB distribution for the
        # speed
        # (i.e. Gaussian distribution with mean=0 and std(=a)=sqrt(3kT/m) for the
        # velocity components
        self.velocities = np.random.normal(0., a, (numParticles,dimension))
        return self.velocities

    ## update functions ##
    def update(self, dT):

        self.updateParticles(dT)
        self.updateVelFor(dT)

        vel2 = self.velocities**2
        vel2sum = np.sum(vel2,axis=None)
        self.energy = (0.5 * m * vel2sum + self.potential) * 10**6 # Energy in Joule
        self.temperature = (2.0/3.0) * (0.5 * m * vel2sum) / (numParticles * kB)  # Paper Verlet 1967

    # update the particle positions
    def updateParticles(self, dT):

        self.positions += self.velocities * dT + 0.5 * (self.force / m) * (dT**2)

        # translate the particles outside of the box
        # +boxSize if positionComponent < 0, -boxSize if positionComponent > 5
        posTranslation = self.positions < 0
        negTranslation = self.positions > boxSize
        self.positions[posTranslation] += boxSize
        self.positions[negTranslation] -= boxSize

    # update the particle velocities and forces
    def updateVelFor(self, dT):

        # TODO calculate forces on particles with the positions (Leo):
        # function:
        # in -> positionVectors (self.position)
        # out -> forceVectors (FORCE)
        #
        # both numParticles by dimension matrices

        self.velocities += 0.5 * (self.force / m) * dT
        self.force, self.potential = Pot.Len_Jones(self.positions)
        self.velocities += 0.5 * (self.force / m) * dT

    # Currently unused!!
    def updateForce(self):

        self.force, self.potential = Pot.Len_Jones(self.positions)



    ## helper functions ##

    # get the particles within a distance imageSize from the box boundaries and translate them to outside
    # the box for force computation to simulate an infinite volume
    #
    # return -> a by d matrix with all the imaged particles
    def getTranslatedImages(self):

        p = self.positions                          # original particle position matrix
        t = np.zeros((np.shape(p)))                 # matrix containing translation values per coordinate
        t[p < imageSize] = boxSize
        t[p > boxSize - imageSize] = -boxSize
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
