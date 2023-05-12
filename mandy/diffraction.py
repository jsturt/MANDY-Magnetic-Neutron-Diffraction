from . import ReadFormFactor as RFF
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
    positions = []
    braggIntensity = []
    braggPosition = []

    # Find MTotal (total magnetic form factor)
    MTotal_list=[]
    # Set unit conversion (p) as unity for now
    p = 1

    # Recover the magnetic supercell
    positions = [row[1:4] for row in crystalObj.pos_df.itertuples()]

    moments = [row[0:4] for row in crystalObj.mom_df.itertuples()]

    # Lambda function to recover the site from a given label
    returnsite = lambda label: [site for site in crystalObj.sites if site.label==label]
    
    # Find the Q = 0 value of the form factor such that it can be normalised
    dictKeys = [site.label for site in crystalObj.sites]
    dictVals = [RFF.form_factor(site.ion_name,0,site.l_s[0],site.l_s[1]) for site in crystalObj.sites ]
    normsDict = dict(zip(dictKeys,dictVals))

    for h, hVal in enumerate(millerIndices[0]):
        for k, kVal in enumerate(millerIndices[1]):
            for l, lVal in enumerate(millerIndices[2]):
                # Progress update
                print('computing ['+str(hVal)+','+str(kVal)+','+str(lVal)+'] ...')
                # Skip the (000) index as it cannot be seen experimentally and leads to div by zero in the calculations.
                if(hVal==0 and lVal==kVal and kVal==hVal):
                    braggIntensity.append( 0.0 )
                    braggPosition.append([hVal,kVal,lVal])
                    continue

                # Clear list for calculation of each miller index
                MTotal_list=[]
                
                # Find q in A^-1 for use with the form factor
                qActual = hVal*crystalObj.b1 + kVal*crystalObj.b2 + lVal*crystalObj.b3

                # Find the magnitude of q, divide by 4pi to match (sin theta / lambda) = (q / 4pi)
                qMag = np.sqrt(np.dot(qActual,qActual)) / (4 * np.pi)
                
                # Calculate the normalised magnetic form factor
                ffList = [RFF.form_factor(returnsite(mom[0])[0].ion_name, qMag, returnsite(mom[0])[0].l_s[0], returnsite(mom[0])[0].l_s[1]) / normsDict[mom[0]] for mom in moments]

                # Apply normalised magnetic form factor
                momentList = [np.array(a[1:4])*b for a,b in zip(moments, ffList)]
                ###momentList = [np.array(a[1:4]) for a,b in zip(moments, ffList)] ### test without form factor

                q = qActual
                #q = [hVal,kVal,lVal]
                                
                # Do vector calculation
                for i,elem in enumerate(positions):
                    MTotalCalc = (p * np.array(momentList[i])) * cmath.exp( complex(0, -np.dot (elem , q) ) ) ### removing 2pi factor
                    # MTotalCalc = (p * np.array(momentList[i])) * cmath.exp( complex(0, -2 * cmath.pi * np.dot (elem , q) ) )
                    MTotal_list.append(MTotalCalc)

                MTotal = sum(MTotal_list)
               
                MTotalOrth = (1/np.dot( q , q )) * np.cross( q , np.cross( MTotal , q ) )

                braggIntensity.append( np.dot( np.conj( MTotalOrth ) , MTotalOrth ) )
                
                braggPosition.append([hVal,kVal,lVal])

    return braggIntensity, braggPosition
                
