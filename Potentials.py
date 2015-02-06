# Calculating potentials and forces
import numpy as np
import variables as var

# Lennard - Jones, sigma and epsilon not yet chosen!
# WARNING: Code is vectorised, i.e. it trades memory for speed. If memory becomes an issue, remove matrixes after they
# have been used, as the code creates multiple size 3*N^2 matrices for calculations

# function Len_Jones.
# input:    positions -> numParticles by dimension matrix holding the potions components of the particles
# output:   force -> numParticles by dimension matrix holding the force components on the particles
#           potentialEnergy -> total potential energy for th particles
def Len_Jones(positions, pImagedParticles = np.empty((0,3))):

    # get the shape of the positions matrix (nP,3) and the imaged particles matrix (nI,3)
    (nP, d) = np.shape(positions)
    (nI, d) = np.shape(pImagedParticles)
    n = nP + nI

    # concatenate the imaged particles behind the positions matrix ie positions -> (n,3)
    positions = np.concatenate((positions, pImagedParticles), axis=0)

    # get the component wise distance (n by 3n) and the distance between particles squared (n by n)
    compDistance = getComponentDistance(positions, n, d)
    r2 = getDistanceSquared(compDistance, n, d)

    # set the distance between imaged particles to infinity to ingnore those interactions
    # for the potential energy
    r2[nP:,nP:] = float("inf")

    # Calculate potential and forces using Leonard-Jones V = 4*eps*((sigma/r)^12-(sigma/r)^6)
    # potentials (forceFactor) is a n by n matrix with on position i,j the potential energy
    # on particle i due to particle j
    potentials = var.eps * ((var.rMin**12 / r2**6) - (var.rMin**6 / r2**3))
    potentialEnergy = 0.5 * np.sum(potentials, axis=None)
    forceFactor = 6.0 * var.eps * ((var.rMin**12 / r2**7) - (var.rMin**6 / r2**4))

    # Create a n x 3*n matrix for the forceFactor
    forceFactor = forceFactor.repeat(d)
    forceFactor = np.reshape(forceFactor,(n,n*d))
    # Calculate the force vectors by multiplying the forceFactor and the difference vectors
    forceMagnitude = forceFactor * compDistance

    # Get the force on each particle by summing contributions of all particles
    force = np.sum(forceMagnitude,axis=0)
    # Reshape to a n x 3 matrix
    force = np.reshape(force, (n,d))

    # remove the imaged particles to get an (np, 3) matrix
    force = force[:nP,:]

    return force, potentialEnergy

# takes an n by 3 matrix with the position components and returns the differences component wise
#
# return -> n by 3n matrix where the (i, 3j + k) element gives the difference of the kth component
# between particle i and j
def getComponentDistance(positions, n, d):

    # Copy the position matrix nt times resulting in a nt by 3nt matrix
    repPos1 = np.kron(np.ones((1,n)),positions)

    # make a length 3N row vector by pasting the position vectors behind each other then paste
    # n times in another n by 3N matrix
    rowPos = np.reshape(positions, (1,n*d))
    repPos2 = np.kron(np.ones((n,1)), rowPos)

    # Subtract the two to get the relative component distances between the coordinates
    return repPos2 - repPos1

# takes an n by 3n matrix with the componentwise differences (see function above)
#
# return -> n by n matrix where the (i,j) element gives the distance squared between particle
# i and j. The (i,j) element is infinity if particle i is particle j (on the diagonal)
def getDistanceSquared(compDistance, n, d):
    # Get the sum of the square distances (with some reshapes to sum) r2 is the
    # distance squared
    compDistanceSquare = np.reshape(compDistance**2,(n**2,d))
    r2 = np.sum(compDistanceSquare, axis=1)
    r2 = np.reshape(r2,(n,n))

    # Replace zeros by infinity to avoid division by zero
    r2[(r2 == 0)] = float("inf")

    return r2

# function constant gravitation field
#
# input:    positions -> numParticles by dimension matrix holding the potions components of the particles
# output:   force -> numParticles by dimension matrix holding the force components on the particles
#           potentialEnergy -> total potential energy for the particles
def Gravity(positions,g):
    (N, d) = np.shape(positions)
    potential = sum(g * positions[:,2])
    force = np.kron(np.ones((N,1)),[0.0,0.0,-g])
    return force, potential