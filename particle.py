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
        phi = _pi * (2 * _random() - 1)

        px = pt * _cos(phi)
        py = pt * _sin(phi)

        #Initialize the four vector.
        self.vec.SetPx(px)
        self.vec.SetPy(py)
        self.vec.SetPz(pz)
        self.vec.SetE(_sqrt(self.vec.P2() + m**2))

        self.isGood = True
        self.veto = False
        if (m == 0):
            self.isFinalState = True
        else:
            self.isFinalState = isFinalState
        self._epsilon = epsilon
        return


    def PTCuts(self, cut):
        if (self.vec.Pt() <= cut):
            self.veto = True
        return

    
    def EtaCuts(self, cut):
        if (abs(self.vec.Eta()) >= cut):
            self.veto = True
        return

    
class Mother(Particle):
    """Particle class capable of decays."""


    def Decay(self, dM1=0, dM2=0, isFinalState1=False,
                       isFinalState2=False):
        """Decay into daughters in the CM frame."""
        assert (not self.isFinalState), 'Can\'t decay final state particles!'
        assert (self.m != 0), 'Can\'t decay massless particles!'
        assert (dM1 + dM2 <= self.m), 'Daughter masses violate CoE!'
                #less than or equal prevent zero-division in rounding errors.
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
        v1.Boost(boost)
        v2.Boost(boost)
        
        if (v1.M2() - dM1**2 > self._epsilon or
            v2.M2() - dM2**2 > self._epsilon):
            self.isGood = False
            
        return (KnownParticle(v1, dM1, isFinalState1, self._epsilon),
                KnownParticle(v2, dM2, isFinalState2, self._epsilon))
        
    
    

class KnownParticle(Mother):
    """Class for a particle with known four-vector."""
    def __init__(self, vec, m, isFinalState, epsilon=1e-7):
        self.vec = vec
        self.m = m
        if(m == 0):
            self.isFinalState = True
        else:    
            self.isFinalState = isFinalState
        self._epsilon = epsilon
        self.veto = False
        self.isGood = True
        
