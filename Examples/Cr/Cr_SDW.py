"""
Demonstration file : Produce a scan along one axis of the transverse SDW in Chromium.
"""

import numpy as np
import mandy as md
from mandy import mandyCrystal as mc
import matplotlib.pyplot as plt

### setup ###
# modulation parameters:
# SDW wavevector
qSDW = np.array([0,0,0.05])
# number of unit cells to tile in each direction to produce the supercell 
n = np.array([2,2,22])
# u and v define the plane in which the modulation takes place. Here we set v=0 to create a SDW 
u = np.array([1,0,0])
v = 0
# miller indices to be investigated
planes = [ [0] , [1] , [round(0 + (0.025*i) , 3 ) for i in range(100)] ]
# define magnetic sites:
# bcc corner sites
cr = mc.site( [1,0,0] , md.requestNIST.find_L_S('Cr1') , ion_name='Cr1' , label='Cr' )
# bcc centre site
cr0 = mc.site( [-1,0,0] , md.requestNIST.find_L_S('Cr1') , ion_name='Cr1' , label='Cr0' )

### crystal creation ###
# load structural and magnetic information
crystal = mc.mandyCrystal('./Chromium.cif', [cr,cr0])
# produce unit cell 
crystal.build()
# visualise unit cell
crystal.plotCrystal()
# create long-range modulation based upon our chosen parameters
crystal.createModulation( qSDW , u, v , 1 , 1, n )
# visualise supercell
crystal.plotCrystal()

### calculation ###
bI, bP = md.diffraction.magnetic_calc( crystal , planes )

### Plotting ###
bN = [str(word[2]) for word in bP]
fig = plt.figure(figsize=(16,8))
ax = plt.axes()
r = np.arange(len(bI))
ax.bar(r,bI)
ax.set_xticks(r,bN,rotation=90,fontsize=11)
fig.tight_layout()
plt.show()
