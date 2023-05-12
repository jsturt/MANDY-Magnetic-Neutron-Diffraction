"""
Demonstration file : Produce a 2D scan of the Bragg peaks in the k=0 plane of a longitudinal SDW in ferrimagnetic NbFe2. 
"""
import numpy as np 
import mandy as md
from datetime import datetime
from mandy import mandyCrystal as mc
import matplotlib.pyplot as plt

### setup ###
out_filename='./data/NbFe2_2DScan_[2,2,11]'
# modulation parameters
qSDW =  np.array([0,0,0.1])
n    =  np.array([2,2,11] )
u    =  np.array([0,0,1]  )
v    =  0
# Nb sites
nb   = mc.site( [0,0,0]      , md.requestNIST.find_L_S('Nb0') , ion_name='Nb0'  , label='Nb'   )
# fe2a site
fe2a = mc.site( [0,0,1]      , md.requestNIST.find_L_S('Fe3') , ion_name='Fe2'  , label='Fe2a' )
# fe6h site
fe6h = mc.site( [0,0,-2.355] , md.requestNIST.find_L_S('Fe3') , ion_name='Fe2'  , label='Fe6h' )
# grid positions
x1 = -0.4
x2 = 3.4
y1 = y2 = 0.0
z1 = -0.4
z2 = 3.4
step = 0.05
nx = int((x2-x1)/step)
nz = int((z2-z1)/step)

### crystal creation ###
crystal = mc.mandyCrystal('./NbFe2_mp-568901_primitive.cif', [nb,fe2a,fe6h])
crystal.build()
crystal.createModulation( qSDW , u, v , 1 , 1, n )

### calculation ###
# construct grid of miller indicies
indices = [[x1 + step*ix for ix in range(nx)] , [y1] , [z1 + step*iz for iz in range(nz)]]
# perform calculation
intensity, positions = md.diffraction.magnetic_calc(crystal, indices)
# reshape and process output such that it may be plotted
z = np.reshape(np.real(np.array(intensity)) , (nx,nz) )

### output ###
np.save(out_filename+'_'+str(datetime.now())[:-7].replace(':','-').replace(' ','_'), z)

### plotting ###
# clamp intensities above and below 3 std deviations in order to better visualise relevant peaks
n_std = 3
vmax = min(z.max()  , round(np.average(z) + (n_std*np.std(z)), 0))
vmin = max(0        , round(np.average(z) - (n_std*np.std(z)), 0))
# figure parameters
plt.rcParams['figure.figsize'] = [5,10]
plt.imshow(z, cmap='turbo', interpolation='kaiser', vmin=vmin, vmax=vmax)
ax = plt.gca()
ax.invert_yaxis()
ax.grid(which='major', alpha=0.5)
ax.grid(which='minor', alpha=0.025)
ax.set_xticks(np.arange(0,np.arange(x1,x2,step).shape[0],10)-0.5, np.round(np.arange(x1,x2,step*10),2), rotation=90)
ax.set_yticks(np.arange(0,np.arange(z1,z2,step).shape[0],4 )-0.5, np.round(np.arange(z1,z2,step*4),2))
ax.set_xticks(np.arange(0,np.arange(x1,x2,step).shape[0],1 )-0.5, np.round(np.arange(x1,x2,step*1 ),2), rotation=90, minor=True)
ax.set_yticks(np.arange(0,np.arange(z1,z2,step).shape[0],1 )-0.5, minor=True)
ax.set_xlabel('a-dir scan')
ax.set_ylabel('c-dir scan')
plt.colorbar(label='intensity (clamped arb. units)')
plt.gcf().tight_layout()
plt.show()
