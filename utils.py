from __future__ import division
from math import pi as _pi, sqrt as _sqrt, atanh as _atanh, atan2 as _atan2


def Eta(pz, p):
    """Calculate pseudorapidity of a particle from pz and p."""
    try:
        eta = _atanh(pz / p)
    except ZeroDivisionError:
        eta = _atanh(pz / (p + 1e-10))
    return eta


def Phi(py, px):
    """Calculate phi of a particle."""
    phi = _atan2(py, px)
    if( phi < 0 ):
        phi += 2 * _pi
    return phi


def CosTheta(vec1, vec2):
    """Return the cosine of the angle between two Lorentz vectors.
    If four vectors are used, they must be of the form (tt, xx, yy, zz)."""
    if( len(vec1) < 3 or
        len(vec2) < 3 ):
        raise RuntimeError('Vector length < 3.')

    elif( len(vec1) != len(vec2) ):
        raise RuntimeError('Vector lengths do not match.')

    elif( len(vec1) == 3 ): #Vectors are of equal length already
        dot = vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]
        dot /= (_sqrt(vec1[0]**2 + vec1[1]**2 + vec1[2]**2) *
                _sqrt(vec2[0]**2 + vec2[1]**2 + vec2[2]**2))
        
     elif( len(vec1) == 4 ):
        dot = vec1[1]*vec2[1] + vec1[2]*vec2[2] + vec1[3]*vec2[3]
        dot /= (_sqrt(vec1[1]**2 + vec1[2]**2 + vec1[3]**2) *
                _sqrt(vec2[1]**2 + vec2[2]**2 + vec2[3]**2))

     else:
         raise RuntimeError(
             'Error calculating dot product. Check vector length')

     return dot
