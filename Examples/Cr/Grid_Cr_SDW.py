"""
Demonstration file : Produce a 2D scan of the Bragg peaks in the h=0 plane of the transverse SDW in Chromium. 
"""
import numpy as np 
import mandy as md
from mandy import mandyCrystal as mc
from datetime import datetime
import matplotlib.pyplot as plt

### setup ###
# output file name
out_filename='~/Documents/python/test_mandy/diffraction/data/Cr_2DScan_[4,4,22]'
# modulation parameters:
# SDW wavevector
qSDW = np.array([0,0,0.05])
# number of unit cells to tile in each direction to produce the supercell
n = np.array([4,4,22])
# u and v define the plane in which the modulation takes place. Here we set v=0 to create a SDW
u = np.array([1,0, 0])
v = 0
# define magnetic sites:
# bcc corner sites
cr  = mc.site( [ 1,0,0] , md.requestNIST.find_L_S('Cr1') , ion_name='Cr1' , label='Cr'  )
# bcc centre site
cr1 = mc.site( [-1,0,0] , md.requestNIST.find_L_S('Cr1') , ion_name='Cr1' , label='Cr0' )
# define grid of positions to be simulated
z1 = -0.15
y1 = -0.5
x1 = x2 =  0.0
z2 = 2.3
y2 = 1.5
step = 0.025
nz = int((z2-z1)/step)
ny = int((y2-y1)/step)

### crystal creation ###
crystal = mc.mandyCrystal('./Chromium.cif', [cr,cr1])
crystal.build()
crystal.createModulation( qSDW , u, v , 1 , 1, n )

### calculation ###
z = np.empty((ny,nz))
for iy in range(ny):
    for iz in range(nz):
            print((iy*nz + iz) / (nz*ny))
            print([ [x1] , [y1 + iy*step] , [z1 + iz*step] ] )
            result = md.diffraction.magnetic_calc(crystal , np.array( [ [x1] , [y1 + iy*step] , [z1 + iz*step] ] ))[0]
            # result = [np.exp((iy))]
            z[iy,iz]=result[0]

### output ###
# concatenate file name with current date+time such that previous simulations are not clobbered.
np.save(out_filename+'_'+str(datetime.now())[:-7].replace(':','-').replace(' ','_'), z)

### plotting ###
# clamp intensities above and below 3 std deviations in order to better visualise relevant peaks
n_std = 3
vmax = min(z.max()  , round(np.average(z) + (n_std*np.std(z)), 0))
vmin = max(0        , round(np.average(z) - (n_std*np.std(z)), 0))
# figure parameters
plt.rcParams['figure.figsize'] = [16, 9]
plt.imshow(z, cmap='turbo', interpolation='kaiser', vmin=vmin, vmax=vmax)
ax = plt.gca()
ax.invert_yaxis()
ax.grid(which='major', alpha=0.5)
ax.grid(which='minor', alpha=0.025)
ax.set_xticks(np.arange(0,np.arange(z1,z2,step).shape[0],10)-0.5, np.round(np.arange(z1,z2,step*10),2), rotation=90)
ax.set_yticks(np.arange(0,np.arange(y1,y2,step).shape[0],10)-0.5, np.round(np.arange(y1,y2,step*10),2))
ax.set_xticks(np.arange(0,np.arange(z1,z2,step).shape[0],1 )-0.5, np.round(np.arange(z1,z2,step*1 ),2), rotation=90, minor=True)
ax.set_yticks(np.arange(0,np.arange(y1,y2,step).shape[0],1 )-0.5, minor=True)
ax.set_xlabel('c-dir scan')
ax.set_ylabel('b-dir scan')
plt.colorbar(label='intensity (clamped arb. units)')
plt.gcf().tight_layout()
plt.show()
