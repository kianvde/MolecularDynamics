import force_loop as fl
import numpy as np
import Potentials as pot
from datetime import datetime

# generate  positions and imaged positions
a = np.random.random((1000,3))
b = np.random.random((50,3))

# get the forces by both methods and keep track of the time
start1 = datetime.now()
(forcePy, potPy) = pot.Len_Jones(a, b)
tPy = datetime.now() - start1

start2 = datetime.now()
(forceFort, potFort) = pot.Len_Jones_Fortran(a, b)
tFort = datetime.now() - start2

# compare the resulting matrices
equality = np.allclose(forcePy, forceFort, rtol=10e-10)

print "------equality methods with relative tolerance 10e-10------"
print equality
print "-----------------------------------------------------------"
print "Ep(fort): ", potFort, ", Ep(python): ", potPy
print "-----------------------------------------------------------"
print "time(fort): ", tFort, ", Ep_py: ", tPy