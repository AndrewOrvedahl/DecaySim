from __future__ import division, print_function
from utils import *
from math import sqrt as _sqrt, pi as _pi, sin as _sin, cos as _cos
from random import random as _random
from vector import VecFour, VecThree
import sys


class Particle(object):
    """Base class for particle objects."""


    def __init__(self, m, pt, pz=0, isFinalState=False, epsilon=1e-7):

        self.vec = VecFour()
        self.m = m
        self.phi = _pi * (2 * _random() - 1)

        px = pt * _cos(self.phi)
        py = pt * _sin(self.phi)

        #Initialize the four vector.
        self.vec.SetPx(px)
        self.vec.SetPy(py)
        self.vec.SetPz(pz)
        self.vec.SetE(_sqrt(self.vec.P2() + m**2))

        self._isGood = True
        self._veto = False
        if (m == 0):
            self._isFinalState = True
        else:
            self._isFinalState = isFinalState
        self._epsilon = epsilon
        return

class Mother(Particle):
    """Particle class capable of decays."""


    def Decay(self, dM1=0, dM2=0, isFinalState1=True,
                       isFinalState2=False):
        """Decay into daughters in the CM frame."""
        if (self._isFinalState):
            return
        assert (dM1 + dM2 <= self.m), 'Daughter masses violate CoE!'
        phi = _pi * (2 * _random() - 1)
        theta = _pi * _random()

        d1rest = VecFour()
        d2rest = VecFour()

        e1rest = (self.m**2 + dM1**2 - dM2**2)/(2*self.m)
        e2rest = self.m - e1rest
        pDaughters = _sqrt(e1rest**2 - dM1**2)

        px = pDaughters * _sin(theta) * _cos(phi)
        py = pDaughters * _sin(theta) * _sin(phi)
        pz = pDaughters * _cos(theta)
        
        v1 = VecFour(px, py, pz, e1rest)
        v2 = VecFour(-px, -py, -pz, e2rest)

        boost = self.vec.BoostVector()
        print(boost)
        v1.Boost(boost)
        v2.Boost(boost)
        
        if (v1.M2() - dM1**2 > self._epsilon or
            v2.M2() - dM2**2 > self._epsilon):
            self._isGood = False
            
        return (KnownParticle(v1, dM1, isFinalState1, self._epsilon),
                KnownParticle(v2, dM2, isFinalState2, self._epsilon))
        

    

class KnownParticle(Mother):
    """Class for a particle with known four-vector."""
    def __init__(self, vec, m, isFinalState, epsilon=1e-7):
        self.vec = vec
        self.m = m
        if(m == 0):
            self._isFinalState = False
        else:    
            self._isFinalState = isFinalState
        self._epsilon = epsilon

