import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import scipy.optimize as sc

import mandy as md


#=======================#===============================================================================================================
#                       #
#   User Parameters:    #
#                       #
#########################

n = np.array([1,1,1])    # number of unit cells used to construct the supercell
qSDW = np.array([0,0,0.1])     # Propagation vector of the magnetic order

# Set up the magnetic site
nameElement = 'Fe2'
# L_element = 0
# S_element = 2.5 

L_element, S_element = md.requestNIST.find_L_S(nameElement)

momentDir = np.array([0,0,1])

# Q = [ [1], [0], [0,1,2,3,4,5,6,7,8,9] ] # Miller indices to be investigated in each direction
Q = [ [0, 1], [0], [-0.1, 0.1, 0.9, 1.1, 1.9, 2.1, 2.9, 3.1] ] # Miller Indicies to be investigated in each direction



#=======================#===============================================================================================================
#                       #
#  Running Simulation:  #
#                       #
#########################

# Build the crystal object
crys = md.mandyCrystal.mandyCrystal(r'NbFe2_mp-568901_primitive.cif')
crys.build()

crys.createSDW(qSDW)

# Perform diffraction calculation
braggIntensity, braggPosition = md.diffraction.magnetic_calc(crys,nameElement,L_element,S_element,Q,momentDir)


#=======================#===============================================================================================================
#                       #
#   Plotting Results:   #
#                       #
#########################

# 3D Plot
#fig = plt.figure(figsize=(12, 12))
#ax = fig.add_subplot(projection='3d')

#def cols(sitename):
#    if(sitename=='Fe2a'):return 'b';
#    if(sitename=='Fe6h'):return 'k';
#    if(sitename=='Nb'):return 'y';

#[ax.scatter(row[1],row[2],row[3],color=cols(row[0]),label='{}'.format(row[0])) for row in crys.pos_df.itertuples()]
#ax.view_init(elev=90, azim=-90)
#scaling = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
#ax.auto_scale_xyz(*[[np.min(scaling), np.max(scaling)]]*3)
#plt.show()

scatter_points = px.scatter_3d(crys.pos_df.reset_index(),x='x',y='y',z='z',color='site_name',size_max=2)
scatter_points.update_layout(scene_aspectmode='manual',
                  scene_aspectratio=dict(x=1, y=1, z=10))
scatter_points.add_trace(go.Cone(
    x = np.array(crys.pos_df.x),
    y = np.array(crys.pos_df.y),
    z = np.array(crys.pos_df.z),
    u = np.array(crys.mom_df.m1),
    v = np.array(crys.mom_df.m2),
    w = np.array(crys.mom_df.m3),
    sizemode = "absolute",
    sizeref = 0.6,
    anchor = "tail"))

scatter_points.show()

# Intensity Plot
braggNames = [str(word) for word in braggPosition]
bw = 0.35
r = np.arange(len(braggIntensity))
r1 = r+(bw)

fig2, ax2 = plt.subplots()
ax2.bar(r1, braggIntensity,zorder=3, width=bw,color=(0.9,0.9,1.0,1),edgecolor='black',hatch='//')

# Plotting Configuration
plt.rcParams["figure.figsize"] = (15,8)
plt.grid(zorder=0)
plt.legend(prop={'size': 13},markerscale=0.1)
plt.xticks(r1,braggNames,fontsize = 15)
#plt.xticks(MomentRatios,newPos,fontsize = 15)
plt.yticks(fontsize = 16)
plt.tight_layout()
plt.subplots_adjust(left=0.1, right=0.95, top=0.9, bottom=0.2)
#plt.xlabel('Ratio of moments 2a:6h', fontsize=16)
plt.xlabel('Miller index', fontsize=16)
#plt.ylabel('Bragg intensity at [1,0,0.1]', fontsize=16)
plt.ylabel('Bragg intensity', fontsize=16)

plt.xticks(rotation = 45)
plt.show()



#=======================================================================================================================================




#Scipy minimising
def fitfunc(ratio, nameObj, L_obj, S_Obj, momObj):
    newQ = [ [1], [0], [0.1] ]

    crystalObj = md.mandyCrystal.mandyCrystal(r'NbFe2_mp-568901_primitive.cif')
    crystalObj.build()

    crystalObj.mom_df.at['Fe6h','m3']=ratio
    
    crystalObj.createSDW(qSDW)
    
    I, P = md.diffraction.magnetic_calc(crystalObj, nameObj, L_obj, S_Obj, newQ, momObj)
    return I[0]


#res = sc.minimize(fitfunc,1,args=(nameElement,L_element,S_element,momentDir))
#print(res)




### Miniming

# def magnetisation(moments,positions,nCells,crys):
#     moment_list = [np.array( moments.loc[ row4[0] ] )[2] * math.cos(2 * math.pi * (np.array(row4[3]) + i)* qSDW ) for i in range(n) for row4 in positions.itertuples() ]
#     return sum(moment_list[0:len(crys.atoms)])/(crys.volume)
# print('magnetisation for current configuration: '+str(magnetisation(mom_df,new_pos,n,crystalData)))


# Loop to find optimal moments to minimise the [1,0,0.1] peak
# MomentRatios = [n/4 for n in range(-1*4,1*4+1)]
# Q = [ [1], [0], [0.1] ]
# newBragg = []
# newPos = ["1 : {}".format(word) for word in MomentRatios]

# def minfunc(ratios, newQ, nameObj, L_obj, S_Obj, momObj):
#     for mom in MomentRatios:
#         crystalObj = md.mandyCrystal.mandyCrystal(r'NbFe2_mp-568901_primitive.cif')
#         crystalObj.build()
#         crystalObj.mom_df.at['Fe6h','m3']=mom
#         crystalObj.createSDW(qSDW)
#         I, P = md.diffraction.magnetic_calc(crystalObj, nameObj, L_obj, S_Obj, newQ, momObj)
#         newBragg.append(I[0])
#     return newBragg
# newBragg = minfunc(MomentRatios,[[1],[0],[0.1]],nameElement,L_element,S_element,momentDir)
# plt.bar(MomentRatios,newBragg,zorder=3,width=bw/2,color=(0.9,0.9,1,1),edgecolor='black',hatch=' //')

# #Scipy minimising
# def fitfunc(ratio,name,L,S,positions,moments,Millerindices,kSpaceLengths,kSpaceAngles,momentOrientation):
#     newQ = [ [1], [0], [0.1] ]
#     moments.at['Fe6h','m3']=ratio
#     I, P= magnetic_calc(name,L,S,positions,moments,newQ,kSpaceLengths,kSpaceAngles,momentOrientation)
#     return I[0]
# res = sc.minimize(fitfunc,1,args=(nameElement,L_element,S_element,new_pos,mom_df,Q,reciprocalLengths,reciprocalAngles,momentDir))
# print(res)

# ## Minimising magnetisation to find antiferromagnetic configuration
# def minfunc(ratio):
#     mom_df.at['Fe6h','m3']=ratio
#     print(ratio)
#     return magnetisation(mom_df,new_pos,n,crystalData)
# res = sc.fsolve(minfunc,1)
# print("optimal ratio of moments for AFM is: "+str(res))
