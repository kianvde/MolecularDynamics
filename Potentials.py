# Calculating potentials and forces

# Lennard - Jones, sigma and epsilon not yet chosen!
def Len_Jones(Pos1, Pos2):
    diff = Pos1 - Pos2
    diff_squaresum = sum(diff**2)
    Potential = 4.0 * ( (1.0 / diff_squaresum**6) - (1.0 / diff_squaresum**3) )
    Force_factor = 4.0 * ( (12.0 / diff_squaresum**7) - (6.0 / diff_squaresum**4) )
    Force1 = Force_factor * diff
    Force2 = -Force1
    return Force1, Force2, Potential

# Gravitational potential along z-axis
def Gravity(Pos,Grav_const):
    Potential = Grav_const * Pos[2]
    Force = [0.0,0.0,-Grav_const]
    return Force, Potential
