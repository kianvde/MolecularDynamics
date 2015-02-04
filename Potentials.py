# Calculating potentials and forces
import numpy as np
import variables as var

# Lennard - Jones, sigma and epsilon not yet chosen!
# WARNING: Code is vectorised, i.e. it trades memory for speed. If memory becomes an issue, remove matrices after they
# have been used, as the code creates multiple size 3*N^2 matrices for calculations
def Len_Jones(Pos_mat):
    shape = np.shape(Pos_mat)
    # Copy the position matrix N times (N = shape[0] = number of particles), result = N x 3*N
    Pos_rep = np.kron(np.ones((1,shape[0])),Pos_mat)
    # Create a 'transpose' of the repeated matrix, result = N x 3*N
    Pos_vec = np.reshape(Pos_mat,(1,shape[0]*shape[1]))
    Pos_repT = np.kron(np.ones((shape[0],1)),Pos_vec)
    # Subtract the two to get the relative distances
    diff_mat = Pos_repT - Pos_rep
    # Get the sum of the square distances (includes some awkward reshapes to facilitate summing)
    diff_square = diff_mat**2
    diff_square = np.reshape(diff_square,(shape[0]**2,shape[1]))
    diff_squaresum = np.sum(diff_square,axis=1)
    diff_squaresum = np.reshape(diff_squaresum,(shape[0],shape[0]))
    # Replace zeros by infinity to avoid division by zero
    diff_squaresum = (diff_squaresum == 0).choose(diff_squaresum,float("inf"))
    # Calculate potential and forces using Lennard-Jones
    Potential_vec = var.epsi * ( (var.r_min**6 / diff_squaresum**6) - (var.r_min**6 / diff_squaresum**3) )
    Potential = 0.5 * np.sum(Potential_vec,axis=None)
    Force_factor = var.epsi * ( (12.0 * var.r_min**7 / diff_squaresum**7) - (6.0 * var.r_min**4 / diff_squaresum**4) )
    # Create a N x 3*N matrix for the force factor
    Force_factor = Force_factor.repeat(shape[1])
    Force_factor = np.reshape(Force_factor,(shape[0],shape[0]*shape[1]))
    # Calculate the force vectors by multiplying the force factor and the difference vectors
    Force_mat = Force_factor * diff_mat
    # Get the force on each particle by summing contributions of all particles
    Force = np.sum(Force_mat,axis=0)
    # Reshape to a N x 3 matrix
    Force = np.reshape(Force,shape)
    return Force, Potential


#Gravitational potential along z-axis
def Gravity(Pos_mat,Grav_const):
    shape = np.shape(Pos_mat)
    Potential = sum(Grav_const * Pos_mat[:,2])
    Force = np.kron(np.ones((shape[0],1)),[0.0,0.0,-Grav_const])
    return Force, Potential
