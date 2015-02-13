import numpy as np
import potentials as pot
from datetime import datetime

i = 1
equality = True
while (equality):
    # generate  positions and imaged positions
    a = np.random.random((50,3))
    b = np.random.random((5,3))

    # get the forces by both methods and keep track of the time
    start1 = datetime.now()
    (forcePy, potPy) = pot.lennardJonesVectorized(a, b)
    tPy = datetime.now() - start1

    start2 = datetime.now()
    (forceFort, potFort) = pot.lennardJones(a, b)
    tFort = datetime.now() - start2

    # compare the resulting matrices
    equalMatrices = np.allclose(forcePy, forceFort, rtol=10e-10)
    equalPotential = np.allclose(potPy, potFort, rtol= 0.01)
    equality = equalMatrices & equalPotential

    print "------equality iteration ", i, "------"
    print equalMatrices, equalPotential
    print "-------------------------------------------------------------"
    print "Ep(fort): ", potFort, ", Ep(python): ", potPy
    print "-------------------------------------------------------------"
    print "time(fort): ", tFort, ", Ep_py: ", tPy
    print "-------------------------------------------------------------"
    i += 1