import math
import pathlib
import numpy as np

"""
Sub-routine containing methods for calculating magnetic form factors.
"""


def form_factor_coeffs(element, order):
    """
    Returns the coefficients in the form factor expansion

    Parameters
    ----------
    element : STRING
        Name of the form factor, e.g. 'Fe2'.
    order : STRING
        Specifies the order of the expansion, e.g. '2'.

    Returns
    -------
    LIST
        List of the coefficients.

    """
    
    # reading file
    filename = pathlib.Path(__file__).parent.joinpath("FormFactorData.dat")
    fFactor = open(filename,"rt")
    
    values = [line for i,line in enumerate(fFactor) if all(elem in line for elem in [element,order]) ]

    nextString = values[0].split()

    return [nextString[m].rpartition("=")[2] for m,val in enumerate(nextString) if "." in val]

def jZero(q2,coeff):
    """
    Returns the j0 term in the form factor dipole approximation

    Parameters
    ----------
    q2 : FLOAT 
        |Q|/4pi.
    coeff : LIST
        Coefficients for the zero'th order term.

    Returns
    -------
    FLOAT
        j0 for position q2.

    """
    coeff = [float(elem) for elem in coeff]
    return coeff[0]*math.exp(-coeff[1]*q2*q2) + coeff[2]*math.exp(-coeff[3]*q2*q2) + coeff[4]*math.exp(-coeff[5]*q2*q2) + coeff[6]

def jTwo(q2,coeff):
    """
    Returns the j2 term in the form factor dipole approximation

    Parameters
    ----------
    q2 : FLOAT 
        |Q|/4pi.
    coeff : LIST
        Coefficients for the zero'th order term.

    Returns
    -------
    FLOAT
        j2 for position q2.

    """
    coeff = [float(elem) for elem in coeff]
    return (q2**2)*(coeff[0]*math.exp(-coeff[1]*q2*q2) + coeff[2]*math.exp(-coeff[3]*q2*q2) + coeff[4]*math.exp(-coeff[5]*q2*q2) + coeff[6])

# Vectorising above functions to avoid repeated calls 

j0_vector = np.vectorize(jZero)
j0_vector.excluded.add(1)

j2_vector = np.vectorize(jTwo,excluded=['coeff'])
j2_vector.excluded.add(1)

def form_factor(elem,q,l,s):
    """
    Calculates the form factor in the dipole approximation in the manner specified by the ILL

    Parameters
    ----------
    elem : STRING
        Name of the form factor, e.g. 'Fe2'.
    q : FLOAT
        |Q|/4pi.
    l : FLOAT
        L Quantum number corresponding to elem.
    s : FLOAT
        S Quantum number corresponding to elem.

    Returns
    -------
    FLOAT
        Form factor.

    """
    return ( (l+2*s)*j0_vector(q,form_factor_coeffs(elem, "0")) + l*j2_vector(q, form_factor_coeffs(elem, "2")) )
    
def form_factor_squared(elem,q,l,s):
    """
    Calculates the squared form factor in the dipole approximation in the manner specified by the ILL

    Parameters
    ----------
    elem : STRING
        Name of the form factor, e.g. 'Fe2'.
    q : FLOAT
        |Q|/4pi.
    l : FLOAT
        L Quantum number corresponding to elem.
    s : FLOAT
        S Quantum number corresponding to elem.

    Returns
    -------
    FLOAT
        Squaredorm factor.
    """
    return ( (l+2*s)*j0_vector(q,form_factor_coeffs(elem, "0")) + l*j2_vector(q, form_factor_coeffs(elem, "2")) )**2    



