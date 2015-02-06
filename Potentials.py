# Calculating potentials and forces
import numpy as np
import variables as var

# Leonard - Joenes, sigma and epsilon not yet chosen!
# WARNING: Code is vectorised, i.e. it trades memory for speed. If memory becomes an issue, remove matrixes after they
# have been used, as the code creates multiple size 3*N^2 matrices for calculations

# function Len_Jones.
# input:    positions -> numParticles by dimension matrix holding the potions components of the particles
# output:   force -> numParticles by dimension matrix holding the force components on the particles
#           potentialEnergy -> total potential energy for th particles
def Len_Jones(positions):

    # get the shape of the positions matrix (N,3)
    (N, d) = np.shape(positions)

    # Copy the position matrix N times resulting in a N by 3N matrix
    repPositions1 = np.kron(np.ones((1,N)),positions)

    # make a length 3N row vector by pasting the position vectors behind eachother then paste
    # N times in another N by 3N matrix
    rowPositions = np.reshape(positions,(1,N*d))
    repPositions2 = np.kron(np.ones((N,1)),rowPositions)

    # Subtract the two to get the relative component distances between the coordinates
    componentDistance = repPositions2 - repPositions1

    # Get the sum of the square distances (with some reshapes to sum) r2 is the
    # distance squared
    componentDistanceSquare = np.reshape(componentDistance**2,(N**2,d))
    r2 = np.sum(componentDistanceSquare,axis=1)
    r2 = np.reshape(r2,(N,N))

    # Replace zeros by infinity to avoid division by zero
    r2[(r2 == 0)] = float("inf")

    # Calculate potential and forces using Leonard-Jones V = 4*eps*((sigma/r)^12-(sigma/r)^6)
    # potentials (forceFactor) is a N by N matrix with on position i,j the potential energy
    # on particle i due to particle j
    potentials = var.eps * ((var.rMin**12 / r2**6) - (var.rMin**6 / r2**3))
    potentialEnergy = 0.5 * np.sum(potentials,axis=None)
    forceFactor = 6.0 * var.eps * ((var.rMin**12 / r2**7) - (var.rMin**6 / r2**4))

    # Create a N x 3*N matrix for the forceFactor
    forceFactor = forceFactor.repeat(d)
    forceFactor = np.reshape(forceFactor,(N,N*d))
    # Calculate the force vectors by multiplying the forceFactor and the difference vectors
    forceMagnitude = forceFactor * componentDistance
    # Get the force on each particle by summing contributions of all particles
    force = np.sum(forceMagnitude,axis=0)
    # Reshape to a N x 3 matrix
    force = np.reshape(force, (N,d))

    return force, potentialEnergy


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