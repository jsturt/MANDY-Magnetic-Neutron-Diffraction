import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import scipy.optimize as sc

import mandy as md
from mandy import mandyCrystal as mc

import sys, os

#=======================#===============================================================================================================
#                       #
#   User Parameters:    #
#                       #
#########################

n = np.array([1,1,1])    # number of unit cells used to construct the supercell
qSDW = np.array([0,0,0.1])     # Propagation vector of the magnetic order

# Set up the magnetic sites

nb = mc.site(moment = [0,0,-0.5000], l_s = md.requestNIST.find_L_S('Nb0'), ion_name = 'Nb0', label = 'Nb')
fe2a = mc.site(moment = [0,0,1.0000], l_s = md.requestNIST.find_L_S('Fe2'), ion_name = 'Fe2', label = 'Fe2a')
fe6h = mc.site(moment = [0,0,-2.8350], l_s = md.requestNIST.find_L_S('Fe2'), ion_name = 'Fe2', label='Fe6h')

# nameElement = 'Fe2'
# (L_element,S_element) = md.requestNIST.find_L_S(nameElement)
# momentDir = np.array([0,0,1])

# Q = [ [1], [0], [0,1,2,3,4,5,6,7,8,9] ] # Miller indices to be investigated in each direction
Q = [ [0, 1], [0], [-0.1, 0.1, 0.9, 1.1, 1.9, 2.1, 2.9, 3.1] ] # Miller Indicies to be investigated in each direction
# Q = [ [0, 1], [0], [round(-0.1 + (0.1*step),2) for step in range(33)] ] # Miller Indicies to be investigated in each direction

#=======================#===============================================================================================================
#                       #
#  Running Simulation:  #
#                       #
#########################

# Build the crystal object
crys = mc.mandyCrystal(r'NbFe2_mp-568901_primitive.cif', [nb,fe2a,fe6h])
crys.build()

# Use new createModulation method
u = np.array([0,0,1])
v = 0
#v = np.array([0,1,0]) 

crys.createModulation(qSDW,u,v)

# Perform diffraction calculation
braggIntensity, braggPosition = md.diffraction.magnetic_calc(crys,Q)

bI, bP = md.diffraction.magnetic_calc(crys,Q)

print(bI)

#=======================#===============================================================================================================
#                       #
#   Plotting Results:   #
#                       #
#########################

#crys.plotCrystal(moment_scale = 0.6, aspect_ratio = dict(x=1, y=1, z=10))
crys.plotCrystal(3)

# Intensity Plot
braggNames = [str(word) for word in braggPosition]
bw = 0.35
r = np.arange(len(braggIntensity))
r1 = r+(bw)

fig3, ax3 = plt.subplots()
r1 = np.arange(len(bI))
ax3.bar(r1, bI,zorder=3, width=bw,color=(0.9,0.9,1.0,1),edgecolor='black',hatch='//')

# Plotting Configuration
plt.rcParams["figure.figsize"] = (15,8)
plt.grid(zorder=0)
plt.legend(prop={'size': 13},markerscale=0.1)
plt.xticks(r1,braggNames,fontsize = 11)
#plt.xticks(MomentRatios,newPos,fontsize = 15)
plt.yticks(fontsize = 11)
plt.tight_layout()
plt.subplots_adjust(left=0.1, right=0.95, top=0.9, bottom=0.2)
#plt.xlabel('Ratio of moments 2a:6h', fontsize=16)
plt.xlabel('Miller index', fontsize=12)
#plt.ylabel('Bragg intensity at [1,0,0.1]', fontsize=16)
plt.ylabel('Bragg intensity', fontsize=12)
plt.xticks(rotation = 90)

plt.show()

#=======================================================================================================================================

#Scipy minimising
def fitfunc(ratio):
    newQ = [ [1], [0], [0.1] ]

    fe6h = mc.site(moment = [0,0,ratio], l_s = md.requestNIST.find_L_S('Fe2'), name = 'Fe2', label='Fe6h')
    crystalObj = md.mandyCrystal.mandyCrystal(r'NbFe2_mp-568901_primitive.cif',[nb,fe2a,fe6h])
    crystalObj.build()

    #crystalObj.mom_df.at['Fe6h','m3']=ratio
    
    crystalObj.createModulation(qSDW,u,v)
    
    I, P = md.diffraction.magnetic_calc_vec(crystalObj, newQ)
    return I[0]

# print("[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[")
# res = sc.minimize(fitfunc,1)
# print(res)


def plotfunc(ratio):
    newQ = [ [1], [0], [0.1] ]

    fe6h = mc.site(moment = [0,0,ratio], l_s = md.requestNIST.find_L_S('Fe2'), name = 'Fe2', label='Fe6h')
    crystalObj = md.mandyCrystal.mandyCrystal(r'NbFe2_mp-568901_primitive.cif',[nb,fe2a,fe6h])

    crystalObj.build()

    #crystalObj.mom_df.at['Fe6h','m3']=ratio
    
    crystalObj.createModulation(qSDW,u,v)
    
    return md.diffraction.magnetic_calc_vec(crystalObj, newQ)[0]

# ratios = [-2.5 + (delta * 0.01) for delta in range(20)]
# rbI = []
# for elem in ratios:
#     rbI.append(plotfunc(elem)[0].real)
# fig4, ax4 = plt.subplots()
# r2 = np.arange(len(rbI))
# ax4.bar(r2, rbI,color=(0.9,0.9,1.0,1),edgecolor='black',hatch='//')
# plt.xticks(r1,braggNames,fontsize = 15)
# plt.show()

