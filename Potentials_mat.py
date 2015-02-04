# Calculating potentials and forces
import numpy as np

# Lennard - Jones, sigma and epsilon not yet chosen!
def Len_Jones(Pos_mat):
    shape = np.shape(Pos_mat)
    Pos_rep = np.kron(np.ones((1,shape[0])),Pos_mat)
    Pos_vec = np.reshape(Pos_mat,(1,shape[0]*shape[1]))
    Pos_repT = np.kron(np.ones((shape[0],1)),Pos_vec)
    diff_mat = Pos_rep - Pos_repT
    # print diff_mat
    diff_square = diff_mat**2
    diff_square = np.reshape(diff_square,(shape[0]**2,shape[1]))
    diff_squaresum = np.sum(diff_square,axis=1)
    diff_squaresum = np.reshape(diff_squaresum,(shape[0],shape[0]))
    diff_squaresum = (diff_squaresum == 0).choose(diff_squaresum,float("inf"))
    Potential_vec = 4.0 * ( (1.0 / diff_squaresum**6) - (1.0 / diff_squaresum**3) )
    Potential = 0.5 * np.sum(Potential_vec,axis=None)
    Force_factor = 4.0 * ( (12.0 / diff_squaresum**7) - (6.0 / diff_squaresum**4) )
    Force_factor = Force_factor.repeat(shape[1])
    Force_factor = np.reshape(Force_factor,(shape[0],shape[0]*shape[1]))
    # print Force_factor
    Force_mat = Force_factor * diff_mat
    # print Force_mat
    Force = np.sum(Force_mat,axis=0)
    Force = np.reshape(Force,(shape[0],shape[1]))
    return Force, Potential

# # Lennard - Jones, sigma and epsilon not yet chosen!
# def Len_Jones(Pos_mat):
#     shape = np.shape(Pos_mat)
#     Pos1 = Pos_mat[0,:].tolist()
#     Pos1_rep = Pos1 * (shape[0] - 1)
#     Pos1_np = np.array(Pos1_rep)
#     Pos1_np = np.reshape(Pos1_np,[(shape[0] - 1),shape[1]])
#     diff_mat = Pos1_np - Pos_mat[1:,:]
#     # print diff_mat
#     diff_square = diff_mat**2
#     diff_squaresum = np.sum(diff_square,axis=1)
#     Potential_vec = 4.0 * ( (1.0 / diff_squaresum**6) - (1.0 / diff_squaresum**3) )
#     Potential = sum(Potential_vec)
#     Force_factor = 4.0 * ( (12.0 / diff_squaresum**7) - (6.0 / diff_squaresum**4) )
#     Force_factor = np.reshape(Force_factor,[(shape[0] - 1),1])
#     # print Force_factor
#     Force_mat = Force_factor * diff_mat
#     # print Force_mat
#     Force = np.sum(Force_mat,axis=0)
#     return Force, Potential



# #Gravitational potential along z-axis
# def Gravity(Pos,Grav_const):
#     Potential = Grav_const * Pos[2]
#     Force = [0.0,0.0,-Grav_const]
#     return Force, Potential
