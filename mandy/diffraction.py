from . import ReadFormFactor as RFF
from . import SelectionRule as sr
import numpy as np
import math,cmath

"""
Sub-routine to perform the diffraction calculation.
"""

def magnetic_calc(n,qSDW,name,L,S,positions,moments,millerIndices,kSpaceLengths,kSpaceAngles,momentOrientation):
    """
    Performs the diffraction calculation.
    
    Parameters
    ----------
    n : INT
        Number of unit cells to be used to construct the supercell.
    qSDW : FLOAT
        Periodicity of the SDW.
    name : STRING
        Site name to be used in the simulated.
    L : FLOAT
        L Quantum number corresponding to the site.
    S : FLOAT
        S Quantum number corresponding to the site.
    positions : PANDAS DATAFRAME
        dataframe containing unit cell positions indexed by site name.
    moments : PANDAS DATAFRAME
        dataframe containing the moments associated with each type of site.
    millerIndices : LIST
        List of Miller Indices to be simulated in each direction.
    kSpaceLengths : NUMPY ARRAY (3,1)
        lattice parameters of the crystal.
    kSpaceAngles : NUMPY ARRAY (3,1)
        lattice angles of the crystal.
    momentOrientation : NUMPY ARRAY (3,1)
        Direction along which the moments are aligned.

    Returns
    -------
    braggIntensity : LIST
        List of the intensities corresponding to each position in braggPosition.
    braggPosition : LIST
        List of simulated Miller Indices.

    """
    momentList = [] 
    sfExp = []
    values = []
    braggIntensity = []
    braggPosition = []
    
    # Create the magnetic supercell
    [values.append(np.array(row3[1:4]) + np.array([0,0,i]) ) for i in range(n) for row3 in positions.itertuples() ]

    # Find the Q = 0 value of the form factor such that it can be normalised
    norm = RFF.form_factor_squared(name,0,L,S) 
    
    # Modulate the moment sizes with a cosine wave.
    # momentList = [np.array( moments.loc[ row4[0] ] )[2] * math.cos(2 * math.pi * (np.array(row4[3]) + i)* qSDW ) for i in range(n) for row4 in positions.itertuples() ]           
    # # mom.append( current_moment * cos(2pi * (cpos + supercell index) * wavevector) )
    
    for h, hVal in enumerate(millerIndices[0]):
        for k, kVal in enumerate(millerIndices[1]):
            for l, lVal in enumerate(millerIndices[2]):
                # Skip the (000) index as it cannot be seen experimentally and leads to div by zero in the calculations.
                if(hVal==0 and lVal==kVal and kVal==hVal):
                    continue
                
                # tempMomentList = [np.array( moments.loc[ row[0] ] ) * sr.selection_rule(momentOrientation,( [ hVal, kVal, lVal ] ),kSpaceAngles,kSpaceLengths) for row in positions.itertuples() ]           
                # mom.append( current_moment * selection rule )
                for i in range(n):
                    for row in positions.itertuples():
                        moment = np.array( moments.loc[ row[0] ] )[2]  # Recover moment from dataframe
                        selRule = sr.selection_rule(momentOrientation,( [ hVal, kVal, lVal ] ),kSpaceAngles,kSpaceLengths) # Find selection rule wrt the current Miller index
                        momentList.append( moment*selRule * math.cos(2 * math.pi * (np.array(row[3]) + i) * qSDW ) ) # Implement SDW with selection rule 
                # mom.append( current_moment * selection rule * cos(2pi * (cpos + supercell index) * wavevector )

                [sfExp.append( (np.dot(loopval, [ hVal, kVal, lVal ] )))  for loopval in values ] # Calculating the argument of the structure factor exponential.

                # Calculate the structure factor, taking into account the moment size
                sf_sdw = 0
                for number in range(len(sfExp)): # Run over the 120 values in momentList and sfExp
                    sf_sdw += momentList[number] * cmath.exp( complex(0, -2 * cmath.pi * sfExp[number] ))
                    
                sfExp = []  # Clear sfExp for use in the next loop
                
                # Find q in A^-1 for use with the form factor
                qActual = kSpaceLengths * np.array([ hVal, kVal, lVal ] )
                # Find the magnitude of q, divide by 4pi to match (sin theta / lambda) = (q / 4pi)
                qMag = np.sqrt(np.dot(qActual,qActual)) / (4 * np.pi)
                
                # Append the Bragg peak to the list and modulate it by the selection rule
                braggIntensity.append(round((RFF.form_factor_squared(name,qMag,L,S) / norm) * abs(sf_sdw / n)**2, 8))

                braggPosition.append(  [ hVal, kVal, lVal ]  )
                
    return braggIntensity, braggPosition

