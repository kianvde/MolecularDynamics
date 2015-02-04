# Calculating potentials and forces
import numpy as np

# Lennard - Jones, sigma and epsilon not yet chosen!
def Len_Jones(Pos1, Pos2):
    diff = np.subtract(Pos1,Pos2)
    diff_squaresum = np.sum(np.square(diff))
    Potential = 4.0 * ( (1 / diff_squaresum**6) - (1 / diff_squaresum**3) )
    Force_factor = 4.0 * ( (12.0 / diff_squaresum**7) - (6.0 / diff_squaresum**4) )
    Force1 = np.multiply(diff,[Force_factor,Force_factor,Force_factor])
    Force2 = np.negative(Force1)
    return Force1, Force2, Potential

# Gravitational potential along z-axis
def Gravity(Pos,Grav_const):
    Potential = Grav_const * Pos[2]
    Force = [0,0,-Grav_const]
    return Force, Potential
