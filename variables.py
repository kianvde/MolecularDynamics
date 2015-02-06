### the variables used in the simulation stored in a single file

## imports
import numpy as np

## constants

# number of particles in the system
numParticlesAxis = 10
numParticles = numParticlesAxis**3

deltaT = 1.             # time step
dimension = 3           # dimensionality of the system
boxSize = 5.            # length of the (cubic) box side
T = 300                 # Temperature (in Kelvin)
a = 0.1                 # Maxwell-Boltzmann standard deviation per component sqrt(3kT/m)

# Lennard-Jones
eps = 1.0               # depth of potential well
rMin = 2.0**(1.0/6.0)   # distance at which potential is minimal

# fraction from the edge of the box translated and used for force calculation
forceCalculationFraction = 0.2                  # fraction
imageSize = forceCalculationFraction*boxSize    # size

## particles class
class Particles(object):


    ## initialization ##
    def __init__(self):

        self.initPositions()
        self.initVelocities()

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
            self.positions[j:j+increment, 0] = posAxis       #For every n particles that are on an axis, set coordinates of those n. Coords are in posAxis.
                                                        #Let's say these are the x coords, then we have x0,x1..xn,x0,x1...xn etc.
            if (j%increment**2)==0:               #Here add the 'z' coordinates after n**2
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

    # initialize the particle velocities
    def initVelocities(self):

        # initiate velocities components according to MB distribution for the
        # speed
        # (i.e. Gaussian distribution with mean=0 and std(=a)=sqrt(3kT/m) for the
        # velocity components
        self.velocities = np.random.normal(0., a, (numParticles,dimension))


    ## update functions ##
    def update(self, dT):

        self.updateParticles(dT)
        self.updateVelocities(dT)

    # update the particle positions
    def updateParticles(self, dT):

        self.positions += self.velocities * dT

        # translate the particles outside of the box
        # +boxSize if positionComponent < 0, -boxSize if positionComponent > 5
        posTranslation = self.positions < 0
        negTranslation = self.positions > boxSize
        self.positions[posTranslation] += boxSize
        self.positions[negTranslation] -= boxSize

    # update the particle velocities
    def updateVelocities(self, dT):


        FORCE = 0.

        self.velocities += FORCE * (dT**2)


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