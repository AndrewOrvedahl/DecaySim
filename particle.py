from __future__ import division
from utils import *
from math import sqrt as _sqrt, log as _log, atanh as _atanh
from math import cos as _cos, sin as _sin, acos as _acos, atan2 as _atan2
from math import pi as _pi
from random import random as _random
import sys


class particle(object):
    """Base class for particle objects."""


    def __init__(self, m, pt, pz=0, phi=-4, isFinalState=False, epsilon=1e-5):

        self.m = m
        self.pT = pt
        self.pZ = pz
        self.phi = _pi * (2 * _random() - 1)
        self.pX = self.pT * _cos(self.phi)
        self.pY = self.pT * _sin(self.phi)
        self.p = _sqrt(self.pT**2 + self.pZ**2)
        try:
            self.eta = _atanh(self.pZ / self.p)
        except ValueError:
            self.eta = _atanh(self.pZ / (self.p + 1e-10))
        self.e = _sqrt(self.m**2 + self.p**2) 
        self.vec4 = (self.e, self.pX, self.pY, self.pZ)
        self._isGood = True
        self._isFinalState = isFinalState
        self._epsilon = epsilon
        

class mother(particle):
    
        
    def GetRestFrameFourVectors(self):
        """Return four vectors of daughters in the rest frame."""
        return (self._rest_vec4_1, self._rest_vec4_2)


    def RestFrameDecay(self, m1=0, m2=0):
        """Decay particle in the CM frame."""

        if( self._isFinalState ):
            return

        self.mD1 = m1
        self.mD2 = m2
        e1_rest = (self.m**2 + m1**2 - m2**2)/(2*self.m)
        e2_rest = self.m - e1_rest
        self.p1_rest = _sqrt(e1_rest**2 - m1**2)
        self.p2_rest = -self.p1_rest

        self._phi = _pi * (2 * _random() - 1)
                
        self._theta = _pi * (2 * _random() - 1)
                
        self._pX1_rest = self.p1_rest * _sin(self._theta) * _cos(self._phi)
        self._pY1_rest = self.p1_rest * _sin(self._theta) * _sin(self._phi)
        self._pZ1_rest = self.p1_rest * _cos(self._theta)
        self._pX2_rest = -self._pX1_rest
        self._pY2_rest = -self._pY1_rest
        self._pZ2_rest = -self._pZ1_rest

        self._rest_vec4_1 = [e1_rest, self._pX1_rest, self._pY1_rest,
                             self._pZ1_rest]
        self._rest_vec4_2 = [e2_rest, self._pX2_rest, self._pY2_rest,
                             self._pZ2_rest]
        
        return
        
        
    def SetDaughtersFourVectors(self):
        """Boost daughters to the lab frame"""
        bx = self.pX / self.e
        by = self.pY / self.e
        bz = self.pZ / self.e
        b2 = bx**2 + by**2 + bz**2
        if( b2 >= 1 ):
            self._isGood = False
            print('Particle v >= c!')
        g = self.e / self.m

        d1prod1 = ((bx * self._rest_vec4_1[1]) + (by * self._rest_vec4_1[2]) +
                   (bz * self._rest_vec4_1[3]))
        d2prod1 = ((bx * self._rest_vec4_2[1]) + (by * self._rest_vec4_2[2]) +
                   (bz * self._rest_vec4_2[3]))

        d1prod2 = g * (g * d1prod1 / (1 + g) + self._rest_vec4_1[0])
        d2prod2 = g * (g * d2prod1 / (1 + g) + self._rest_vec4_2[0])

        self._vec4_1 = (g * (self._rest_vec4_1[0] + d1prod1),
                        self._rest_vec4_1[1] + d1prod2 * bx,
                        self._rest_vec4_1[2] + d1prod2 * by,
                        self._rest_vec4_1[3] + d1prod2 * bz)

        self._vec4_2 = (g * (self._rest_vec4_2[0] + d2prod1),
                        self._rest_vec4_2[1] + d2prod2 * bx,
                        self._rest_vec4_2[2] + d2prod2 * by,
                        self._rest_vec4_2[3] + d2prod2 * bz)
        
        self.p1 = _sqrt(self._vec4_1[1]**2 +
                       self._vec4_1[2]**2 +
                       self._vec4_1[3]**2)

        self.p2 = _sqrt(self._vec4_2[1]**2 +
                       self._vec4_2[2]**2 +
                       self._vec4_2[3]**2)

        self.e1 = self._vec4_1[0]
        self.pX1 = self._vec4_1[1]
        self.pY1 = self._vec4_1[2]
        self.pZ1 = self._vec4_1[3]

        self.e2 = self._vec4_2[0]
        self.pX2 = self._vec4_2[1]
        self.pY2 = self._vec4_2[2]
        self.pZ2 = self._vec4_2[3]
        
        if( self.e1**2 - self.p1**2 - self.mD1**2 > self._epsilon or
            self.e2**2 - self.p2**2 - self.mD2**2 > self._epsilon):
            self._isGood = False
            
        return

        
    def GetDaughtersFourVectors(self):
        """Return daughters' four vectors in the lab frame."""
        return (self._vec4_1, self._vec4_2)
    

    def CalculateDaughterObservables(self):
        """Calculate daughter observables."""

        self.phi1 = _atan2(self.pY1, self.pX1)
        self.phi2 = _atan2(self.pY2, self.pX2)
        self.dPhi = self.phi1 - self.phi2
        if( self.dPhi >= _pi ):
            self.dPhi -= 2 * _pi
        if( self.dPhi < -_pi ):
            self.dPhi += 2 * _pi
        self.pT1 = _sqrt(self.pX1**2 + self.pY1**2)
        self.pT2 = _sqrt(self.pX2**2 + self.pY2**2)

        try:
            self.eta1 = _atanh(self.pZ1 / self.p1)
        except ValueError:
            self.eta1 = _atanh(self.pZ1 / (self.p1 + 1e-10))
        try:
            self.eta2 = _atanh(self.pZ2 / self.p2)
        except ValueError:
            self.eta2 = _atanh(self.pZ2 / (self.p2 + 1e-10))
        self.dEta = abs(self.eta1 - self.eta2)

        self.dR = _sqrt(self.dEta**2 + self.dPhi**2)
        self.cosTheta = ((self.pX1 * self.pX2 +
                          self.pY1 * self.pY2 +
                          self.pZ1 * self.pZ2) /
                         (self.p1 * self.p2))


        if( abs(self.pX1 + self.pX2 - self.pX) > self._epsilon or
            abs(self.pY1 + self.pY2 - self.pY) > self._epsilon or
            abs(self.pZ1 + self.pZ2 - self.pZ) > self._epsilon ):
            self._isGood = False
        return


    def DaughterPTCuts(self, pT):
        """Cut events that don't meet pT threshold.""" 
        if( self.pT1 < pT or self.pT2 < pT ):
            self._isGood = False


    def DaughterEtaCuts(self, eta):
        """Cut events that don't meet pseudorapidity threshold.""" 
        if( self.eta1 > eta or self.eta2 > eta ):
            self._isGood = False


    def MotherPTCuts(self, pT):
        """Cut events that don't meet pT threshold."""
        if( self.pT < pT ):
            self._isGood = False


    def MotherEtaCuts(self, eta):
        """Cut events that don't meet pseudorapidity threshold."""
        if( self.eta > eta ):
            self._isGood = False
            

    def GetDaughters(self, isFinalState1, isFinalState2):
        """Return a tuple of the daughter particles."""
        return(particle(self.mD1, self.pT1, pz=self._vec4_1[3], phi=self.phi1,
                        isFinalState=isFinalState1),
               particle(self.mD2, self.PT2, pz=self._vec4_2[3], phi=self.phi2,
                        isFinalState=isFinalState2))


    def GoodEvent(self):
        """Check whether event is good."""
        return self._isGood


    def __str__(self):
        """Return a string in pseudo-LHE format."""
        string = '<evt>\npX\t pY\t pZ\t E\t M\n'
        string += ('%7.4f %7.4f %7.4f %7.4f %7.4f\n' %
                   (self.pX, self.pY, self.pZ, self.e, self.m))
        try:
            string += ('%7.4f %7.4f %7.4f %7.4f %7.4f\n' %
                       (self.pX1, self.pY1, self.pZ1, self.e1, self.mD1))
            string += ('%7.4f %7.4f %7.4f %7.4f %7.4f\n' %
                       (self.pX2, self.pY2, self.pZ2, self.e2, self.mD2))
        except AttributeError:
            pass
        string += '<\evt>'
        
        return string
