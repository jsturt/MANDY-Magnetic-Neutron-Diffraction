import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.express as px
import math, cmath
import pandas as pd
"""
# Arguments:
	u: vector
	v: vector
	m1: scalar
	m2: scalar
	Rl: array of vectors
	rj: array of vectors
	k: vector
"""

# Define the number of unit cells in use
n = [1,11,1]

# Define propagation vector
k = np.array([0,0.1,0])

# Define Miller Index to be investigated
Q = np.array([1,0,1])

# Define the plane of modulation for a SDW
# u = np.array([0,0,1])
# v = np.array([0,0,0])
# Define the plane of modulation for a helical structure
u = np.array([0,0,1])
v = np.array([1j,0,0])
#Define the plane of modulation for a cycloidal structure
# u = np.array([0,1,0])
# v = np.array([0,0,1j])


# Define moment components in the u and v directions 
m1 = 1
m2 = 1

# Define positions for a simple cubic lattice with one atom per unit cell
Rl = [np.array([i,j,k]) for i in range(n[0]) for j in range(n[1]) for k in range(n[2])]
rj = [np.array([0,0,0]),np.array([0.5,0.5,0.5])]


#------------------------------------------------------------------------------------------------------
# MODULATION


# moments = []
# for R in Rl:
# 	for r in rj:
# 		moment = (m1*u + 1j*m2*v) * cmath.exp( complex(0, -2*cmath.pi * (np.dot(k,R) + np.dot(k,r)) ) )
# 		moments.append(moment)


moments = []
for R in Rl:
    for r in rj:
        moment = [(m1*u + 1j*m2*v) * cmath.exp( complex(0, -2*cmath.pi * (np.dot(vec,R) + np.dot(vec,r)) ) ) for vec in [k,-1*k]]
        moments.append(0.5*sum(moment))
        print('complex test : ' + str(any(np.iscomplex(0.5*sum(moment)))))

# no i version
# moments = []
# for R in Rl:
#     for r in rj:
#         moment = [(m1*u - m2*v) * cmath.exp( complex(0, -2*cmath.pi * (np.dot(vec,R) + np.dot(vec,r)) ) ) for vec in [k,-1*k]]
#         moments.append(0.5*sum(moment))
#         print('complex test : ' + str(any(np.iscomplex(0.5*sum(moment)))))



#------------------------------------------------------------------------------------------------------
# INTENSITY

def intensity(h):
    # For now set form factor & the unit coversion (p) as unity
    formfactor = 1
    p = 1

    # Find MTotal (total magnetic form factor)
    MTotal_list = []
    index = 0
    for R in Rl:
        for r in rj:
            moment = (p * formfactor * moments[index]) * cmath.exp( complex(0, -2*cmath.pi * (np.dot(h,R) + np.dot(h,r)) ) )
            MTotal_list.append(moment)
            index += 1

    MTotal = sum(MTotal_list)
    # print('MTotal:')
    # print(MTotal)

    # Find Mtotal_orth
    MTotalOrth = (1/np.dot(h,h)) * np.cross(h,np.cross(MTotal,h))
    # print('MTotalOrth')
    # print(MTotalOrth)

    # Find intensity
    return np.dot(np.conj(MTotalOrth), MTotalOrth)

#print(intensity(Q))


#------------------------------------------------------------------------------------------------------
# LOOPS :- Testing intensity behaviours

braggIntensity = []
braggLabels = []
xaxis = []

for i in range(30):
    xaxis.append(0.1*(i+1))
    braggLabels.append(str(np.array([1,0,0.1*(i+1)])))
    braggIntensity.append(np.real(intensity(np.array([1,0,0.1*(i+1)]))))


#------------------------------------------------------------------------------------------------------
# PLOTTING, doesn't matter if this code is crap.

pos = [R+r for R in Rl for r in rj]


print('Test if pos & Rl are same')
print(all([(Rl[ind]==pos[ind]).all() for ind in range(len(Rl))]))

df = pd.DataFrame(pos, columns=['x','y','z'])
mdf = pd.DataFrame(moments, columns=['x','y','z'])

def plotCrystal(moment_scale = 0.2, aspect_ratio = dict(x=1, y=1, z=1)):
    scatter_points = px.scatter_3d(df,x='x',y='y',z='z',size_max=2)
    scatter_points.update_layout(
                        scene_aspectmode='manual',
                        scene_aspectratio=aspect_ratio)
    scatter_points.add_trace(go.Cone(
        x = np.array(df.x),
        y = np.array(df.y),
        z = np.array(df.z),
        u = np.real(mdf.x),
        v = np.real(mdf.y),
        w = np.real(mdf.z),
        sizemode = "absolute",
        sizeref = moment_scale,
        anchor = "tail",
        colorscale = "blues",
        hoverinfo = 'skip'))
    scatter_points.update_traces(showscale=False, selector=dict(type='cone'))
    scatter_points.show()

plotCrystal()

fig2, ax2 = plt.subplots(figsize=(20, 8))
ax2.bar(xaxis, braggIntensity,zorder=3, width=0.05,color=(0.9,0.9,1.0,1),edgecolor='black',hatch='//')
plt.xticks(xaxis,braggLabels,fontsize = 15, rotation=45) 
plt.tight_layout()
plt.show()