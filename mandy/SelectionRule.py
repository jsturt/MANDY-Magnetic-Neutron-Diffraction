import numpy as np

"""
Sub-routine to enable the selection rule to be applied to intensities.
"""

def theta(one,two):
    """
    Method to returnt the angle between two 3D vectors

    Parameters
    ----------
    one : NUMPY ARRAY (3,1)
        First vector.
    two : NUMPY ARRAY (3,1)
        Second vector.

    Returns
    -------
    FLOAT
        Angle between vector one and two.

    """
    return np.arccos( np.dot(one,two) / ( np.sqrt(np.dot(one,one)) * np.sqrt(np.dot(two,two)) ) )

def selection_rule(inA,inB,ang,sid):
    """
    Method for implementing the selection rule via rotation matrices.

    Parameters
    ----------
    inA : NUMPY ARRAY (3,1)
        The first lattice vector to be considered.
    inB : NUMPY ARRAY (3,1)
        The second lattice vector to be considered.
    ang : NUMPY ARRAY (3,1)
        Lattice angle parameters.
    sid : NUMPY ARRAY (3,1)
        Lattice length parameters.
                            
    Returns
    -------
    FLOAT
        Sin squared value of the true angle between the two vectors.

    """
    threshold=0.005
    if(np.sum(inA)<threshold or np.sum(inB)<threshold):return 0     # Prevent div by zero for small vectors entered.
    # Convert angles to radians
    ang = np.radians(ang)
    # Rotation and scaling of the axes
    a = np.array([0,sid[0],0])
    b = sid[1] * np.array( [np.sin(ang[2])*np.sin(ang[0]) , np.cos(ang[2])*np.sin(ang[0]) , np.cos(ang[0])] )
    c = sid[2] * np.array([0,np.cos(ang[1]),np.sin(ang[1])])
    # Calculate the vectors in cartesian basis
    A = inA[0]*a + inA[1]*b + inA[2]*c
    B = inB[0]*a + inB[1]*b + inB[2]*c
    t = theta(A,B)
    return np.sin(t)**2

# Alternative method tried to find the selection rule.
# def sel_rule2(inA,inB,ang,sid):
#     ang = np.radians(ang)
#     alpha = np.arccos( ( np.cos(ang[1])*np.cos(ang[2]) - np.cos(ang[0]) ) / ( np.sin(ang[1])*np.sin(ang[2]) ) )
#     A = [ sid[0] , sid[1]*np.cos(ang[2]) , sid[2]*np.cos(ang[1]) ]
#     B = [ 0 , sid[1]*np.sin(ang[2]) , - sid[2]*np.sin(ang[1])*np.cos(alpha) ]
#     C = [ 0 , 0 , sid[2]*np.sin(ang[1])*np.sin(alpha) ]
#     M = np.array((A,B,C))
#     outA = np.dot(M,inA)
#     outB = np.dot(M,inB)
#     return np.sin(theta(outA,outB))**2
    

#print("Using the given lattice: theta =  ", theta(A,B)*360/(2*np.pi), " degrees")
#print("In a square lattice: theta = ",theta(inA,inB)*360/(2*np.pi), " degrees")
#
#print("Selection rule value: ",selection_rule(inA,inB,ang,sid))