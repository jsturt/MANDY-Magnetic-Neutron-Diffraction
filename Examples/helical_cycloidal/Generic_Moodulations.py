"""
Demonstration file: Produce cycloidal and helical magnetic structures via the 'rotate' parameter.
"""

import numpy as np 
import mandy as md
from mandy import mandyCrystal as mc
import matplotlib.pyplot as plt

### setup ###

# modulation parameters:
# modulation wavevector.
qMOD = np.array([0,0,0.2])
# number of unit cells to tile in each direction to produce the supercell.
n = np.array([1,1,10])
# u and v define the plane in which the modulation takes place. 
# 	- Here we use non-zero v to create modulation in more than one direction (as would be the case in a SDW).
#	- Choose u,v in plane with the wavevector to produce a cycloidal modulation.
u_cycloidal = np.array([0,0,1])
v_cycloidal = np.array([0,1,0])
# define magnetic site:
cr = mc.site(moment = [0,0,1] , l_s = md.requestNIST.find_L_S('Cr0') , ion_name = 'Cr' , label = 'Cr' )

### crystal creation ###

# load structural and magnetic information.
crystal = mc.mandyCrystal('Chromium.cif',[cr])
# produce unit cell.
crystal.build()
# create long-range modulation based upon our chosen parameters.
# 	- We set rotate=True. This means that only the magnitude of our moments [cr] are taken into account, 
#		the moments are rotated to match the modulation instead of being multiplied with it.
crystal.createModulation(qMOD,u_cycloidal,v_cycloidal,1,1,n,rotate=True)
# visualise supercell.
crystal.plotCrystal(moment_scale=4)

# reset crystal object to plot the next example.
crystal.build()
# choose u,v orthogonal to the wavevector to produce a helical modulation.
u_helical = np.array([1,0,0])
v_helical = np.array([0,1,0])
crystal.createModulation(qMOD,u_helical,v_helical,1,1,n,rotate=True)
# visualise supercell.
crystal.plotCrystal(moment_scale=4)

# Diffraction calculations can then be performed as usual.