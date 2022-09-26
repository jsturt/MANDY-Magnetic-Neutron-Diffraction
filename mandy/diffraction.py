from . import ReadFormFactor as RFF
from . import SelectionRule as sr
import numpy as np
import math,cmath

"""
Sub-routine to perform the diffraction calculation.
"""

def magnetic_calc(crystalObj,millerIndices):
    """
    Performs the diffraction calculation.
    
    Parameters
    ----------
    crystalObj : mandyCrystal object
        Object which contains the crystal on which the simulation is conducted.
    millerIndices : LIST
        List of Miller Indices to be simulated in each direction.

    Returns
    -------
    braggIntensity : LIST
        List of the intensities corresponding to each position in braggPosition.
    braggPosition : LIST
        List of simulated Miller Indices.

    """
    momentList = []
    values = []
    braggIntensity = []
    braggPosition = []

    # Recover the magnetic supercell
    values = [row[1:4] for row in crystalObj.pos_df.itertuples()]

    moments = [row[0:4] for row in crystalObj.mom_df.itertuples()]

    # Lambda function to recover the site from a given label
    returnsite = lambda label: [site for site in crystalObj.sites if site.label==label]
    
    # Find the Q = 0 value of the form factor such that it can be normalised
    dictKeys = [site.label for site in crystalObj.sites]
    dictVals = [RFF.form_factor(site.name,0,site.l_s[0],site.l_s[1]) for site in crystalObj.sites ]
    normsDict = dict(zip(dictKeys,dictVals))

    for h, hVal in enumerate(millerIndices[0]):
        for k, kVal in enumerate(millerIndices[1]):
            for l, lVal in enumerate(millerIndices[2]):
                # Skip the (000) index as it cannot be seen experimentally and leads to div by zero in the calculations.
                if(hVal==0 and lVal==kVal and kVal==hVal):
                    continue

                # Find q in A^-1 for use with the form factor
                qActual = hVal*crystalObj.b1 + kVal*crystalObj.b2 + lVal*crystalObj.b3

                # Find the magnitude of q, divide by 4pi to match (sin theta / lambda) = (q / 4pi)
                qMag = np.sqrt(np.dot(qActual,qActual)) / (4 * np.pi)

                # Find selection rule at the current Miller index
                selRule = [sr.selection_rule(mom[1:4],( [ hVal, kVal, lVal ] ),crystalObj.reciprocalAngles,crystalObj.reciprocalLengths) for mom in moments]
                # Apply the selection rule
                momentList = [np.array(a[1:4]) * b for a,b in zip(moments,selRule)] 
                # Calculate the normalised magnetic form factor
                ffList = [RFF.form_factor(returnsite(mom[0])[0].name, qMag, returnsite(mom[0])[0].l_s[0], returnsite(mom[0])[0].l_s[1]) / normsDict[mom[0]] for mom in moments]
                # Apply normalised magnetic form factor
                momentList = [a*b for a,b in zip(momentList, ffList)]

                # Calculating the argument of the structure factor exponential.
                sfExp = [(np.dot(loopval, [ hVal, kVal, lVal ] ))  for loopval in values ]
                
                # Calculate the structure factor, taking into account the moment size
                sf_sdw = 0
                for number in range(len(sfExp)): # Run over the 120 values in momentList and sfExp
                    sf_sdw += momentList[number] * cmath.exp( complex(0, -2 * cmath.pi * sfExp[number] ))

                sfExp = []  # Clear sfExp for use in the next loop
                                
                # Append the Bragg peak to the list   
                braggIntensity.append( np.real( (sf_sdw).dot(np.conj(sf_sdw) ) ) )

                braggPosition.append(  [ hVal, kVal, lVal ]  )
                
    return braggIntensity, braggPosition

