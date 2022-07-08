from . import ReadFormFactor as RFF
from . import SelectionRule as sr
import numpy as np
import math,cmath

"""
Sub-routine to perform the diffraction calculation.
"""

def magnetic_calc(crystalObj,name,L,S,millerIndices,momentOrientation):
    """
    Performs the diffraction calculation.
    
    Parameters
    ----------
    name : STRING
        Site name to be used in the simulated.
    L : FLOAT
        L Quantum number corresponding to the site.
    S : FLOAT
        S Quantum number corresponding to the site.
    millerIndices : LIST
        List of Miller Indices to be simulated in each direction.
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
    
    # Recover the magnetic supercell
    values = [row[1:4] for row in crystalObj.pos_df.itertuples()]

    moments = [row[1:4] for row in crystalObj.mom_df.itertuples()]

    # Find the Q = 0 value of the form factor such that it can be normalised
    norm = RFF.form_factor_squared(name,0,L,S) 

    
    for h, hVal in enumerate(millerIndices[0]):
        for k, kVal in enumerate(millerIndices[1]):
            for l, lVal in enumerate(millerIndices[2]):
                # Skip the (000) index as it cannot be seen experimentally and leads to div by zero in the calculations.
                if(hVal==0 and lVal==kVal and kVal==hVal):
                    continue

                                   
                selRule = sr.selection_rule(momentOrientation,( [ hVal, kVal, lVal ] ),crystalObj.reciprocalAngles,crystalObj.reciprocalLengths) # Find selection rule wrt the current Miller index
                
                momentList = [np.array(mom) * selRule for mom in moments]
                
                [sfExp.append( (np.dot(loopval, [ hVal, kVal, lVal ] )))  for loopval in values ] # Calculating the argument of the structure factor exponential.

                # Calculate the structure factor, taking into account the moment size
                sf_sdw = 0
                for number in range(len(sfExp)): # Run over the 120 values in momentList and sfExp
                    sf_sdw += momentList[number] * cmath.exp( complex(0, -2 * cmath.pi * sfExp[number] ))
                
                sfExp = []  # Clear sfExp for use in the next loop
                
                # Find q in A^-1 for use with the form factor
                qActual = crystalObj.reciprocalLengths * np.array([ hVal, kVal, lVal ] )
                # Find the magnitude of q, divide by 4pi to match (sin theta / lambda) = (q / 4pi)
                qMag = np.sqrt(np.dot(qActual,qActual)) / (4 * np.pi)
                
                # Append the Bragg peak to the list and modulate it by the selection rule
                braggIntensity.append( (RFF.form_factor_squared(name,qMag,L,S) / norm) * np.dot(sf_sdw / crystalObj.n, sf_sdw / crystalObj.n) )

                braggPosition.append(  [ hVal, kVal, lVal ]  )
                
    return braggIntensity, braggPosition

