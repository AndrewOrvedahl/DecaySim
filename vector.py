from __future__ import division, print_function
from math import sqrt as _sqrt, atan2 as _atan2, log as _log, pi as _pi

#Force compatibility with python 2 and 3.
try:
    xrange
except NameError:
    xrange = range

class VecThree(object):
    """Simple three-vector implementation. Every method for retrieving or
    setting a value is implemented for both momentum and position. The 
    methods are interchangeable."""
    def __init__(self, x=0, y=0, z=0):
        self._x = x
        self._y = y
        self._z = z
        return

    def P2(self):
        """Return the squared magnitude of momentum."""
        return self._x**2 + self._y**2 + self._z**2


    def P(self):
        """Return the magnitude of the momentum."""
        return _sqrt(self._x**2 + self._y**2 + self._z**2)


    def Pt(self):
        """Return the transverse momentum of the particle."""
        return _sqrt(self._x**2 + self._y**2)


    def X(self):
        """Return the x-component."""
        return self._x


    def Px(self):
        """Return the px-component."""
        return self._x


    def Y(self):
        return self._y
    

    def Py(self):
        """Return the py-component."""
        return self._y


    def Z(self):
        """Return the z-component."""
        return self._z


    def Pz(self):
        """Return the pz-component."""
        return self._z


    def SetX(self, x):
        """Set the x-component."""
        self._x = x


    def SetPx(self, x):
        """Set the x-component."""
        self._x = x


    def SetY(self, y):
        """Set the y-component."""
        self._y = y


    def SetPy(self, y):
        """Set the y-component."""
        self._y = y


    def SetZ(self, z):
        """Set the z-component."""
        self._z = z


    def SetPz(self, z):
        """Set the z-component."""
        self._z = z


    def Phi(self):
        """Return Phi of the particle."""
        return _atan2(self._y, self._x)


    def Eta(self):
        """Return pseudorapidity of the particle."""
        p = _sqrt(self._x**2 + self._y**2 + self._z**2)
        if (p == 0):
            return 0.0
        if (p == self._z and p > 0):
            return 1e72
        if (p == self._z and p < 0):
            return -1e72
        return 0.5*_log((p + self._z) / (p - self._z))


    def DeltaPhi(self, v):
        """Return delta phi of it and some other particle."""
        dPhi = self.Phi() - v.Phi()
        if (dPhi >= _pi):
            dPhi -= 2 * _pi
        if (dPhi < -_pi):
            dPhi += 2 * _pi
        return dPhi
    
    
    def DeltaEta(self, v):
        """Return delta eta of it and some other vector."""
        return abs(self.Eta() - v.Eta())


    def DeltaR(self, v):
        """Return delta eta of it and some other vector."""
        return _sqrt(self.DeltaEta(v)**2 + self.DeltaPhi(v)**2)


    def Dot3(self, v):
        """Return the dot product of it's and another vector's 
        positions/momenta."""
        return (self.X() * v.X() + self.Y() * v.Y() + self.Z() * v.Z())


    def Unit3(self):
        """Return a unit three-vector."""
        mag = self.P()
        return VecThree(self._x / mag,
                        self._y / mag,
                        self._z / mag)


    def CosTheta(self, v):
        """Return the cosine of the angle between two vectors."""
        return self.Unit3().Dot3(v.Unit3())


    def __add__(self, v):
        """Return the vector sum."""
        return VecThree(self.X() + v.X(), self.Y() + v.Y(), self.Z() + v.Z())


    def __sub__(self, v):
        """Return the vector difference."""
        return VecThree(self.X() - v.X(), self.Y() - v.Y(), self.Z() - v.Z())


    def __mul__(self, v):
        """Multiply self by a scalar. This should not be used with a VecFour."""
        return VecThree(self.X()*v, self.Y()*v, self.Z()*v)


    def __truediv__(self, v):
        """Divide self by a scalar. This should mot be used with a VecFour."""
        return VecThree(self.X()/v, self.Y()/v, self.Z()/v)


    def __eq__(self, v):
        """Check equality of two vectors."""
        if (self.X() == v.X() and
            self.Y() == v.Y() and
            self.Z() == v.Z()):
            return True
        return False


    def __getitem__(self, key):
        """Allows acces to values via indexing."""
        if (key == 0 or key == 'x' or key == 'px'):
            return self.X()
        if (key == 1 or key == 'y' or key == 'py'):
            return self.Y()
        if (key == 2 or key == 'z' or key == 'pz'):
            return self.Z()
        else:
            raise AttributeError


    def __repr__(self):
        return '(%.4f, %.4f, %.4f)' % (self.X(), self.Y(), self.Z())


    def __str__(self):
        """Represents vector as a string."""
        return '(%.4f, %.4f, %.4f)' % (self.X(), self.Y(), self.Z())

    
    def __imul__(self, scalar):
        """Return product of self and a scalar. This should not be used with
        a VecFour."""
        self.SetX(self.X() * scalar)
        self.SetY(self.Y() * scalar)
        self.SetZ(self.Z() * scalar)
        return self


    def __iadd__(self, v):
        """Adds another vector to self."""
        self.SetX(self.X() + v.X())
        self.SetY(self.Y() + v.Y())
        self.SetZ(self.Z() + v.Z())
        

    def __neg__(self):
        """Flip sign of spcial/momentum components."""
        self.SetX(self.X() * -1.)
        self.SetY(self.Y() * -1.)
        self.SetZ(self.Z() * -1.)


    def Generator(self):
        for i in xrange(3):
            yield self[i]
        
class VecFour(VecThree):
    """Simple four-vector implementation. Every method for retrieving or
    setting a value is implemented for both momentum and position. The 
    methods are interchangeable."""
    def __init__(self, x=0, y=0, z=0, t=0):
        super(VecFour, self).__init__(x, y, z)
        self._t = t
        return


    def T(self):
        """Return the t-component."""
        return self._t


    def E(self):
        """Return the e-component."""
        return self._t


    def M2(self):
        """Return E^2 - P^2."""
        return self.T()**2 - (self.X()**2 + self.Y()**2 + self.Z()**2)
    

    def SetT(self, t):
        """Set the t-component."""
        self._t = t


    def SetE(self, e):
        """Set the e-component."""
        self._t = e


    def Dot4(self, v):
        """Return E1*E2 - P1*P2."""
        return self.T() * v.T() - self.Dot3(v)


    def Unit4(self):
        """Return a unit four-vector."""
        mag = _sqrt(self._x**2 + self._y**2 + self._z**2 + self._t**2)
        return VecFour(self._x / mag,
                       self._y / mag,
                       self._z / mag,
                       self._t / mag)

    
    def BoostVector(self):
        """Return a three-vector conaining Beta-x, Beta-y, and Beta-z."""
        bx = self._x / self._t
        by = self._y / self._t
        bz = self._z / self._t
        return VecThree(bx, by, bz)


    def Boost(self, v):
        """Boost the vector to the rest-frame of some other vector. It accepts
        a VecThree of the components of beta."""
        bx = v.X()
        by = v.Y()
        bz = v.Z()
        b2 = bx**2 + by**2 + bz**2
        g = 1. / _sqrt(1. - b2)
        bp = self.Dot3(v)

        if (b2 > 0):
            g2 = (g - 1.)/b2
        else:
            g2 = 0.

        self.SetX(self.X() + g2*bp*bx + g*bx*self.T())
        self.SetY(self.Y() + g2*bp*by + g*by*self.T())
        self.SetZ(self.Z() + g2*bp*bz + g*bz*self.T())
        self.SetT(g*(self.T() + bp))
        return


    def __add__(self, v):
        """Return the vector sum."""
        return VecFour(self.X() + v.X(), self.Y() + v.Y(),
                       self.Z() + v.Z(), self.T() + v.T())


    def __sub__(self, v):
        """Return the difference of the vectors' components."""
        return VecThree(self.X() - v.X(), self.Y() - v.Y(), self.Z() - v.Z())


    def __eq__(self, v):
        """Check equality of two vectors."""
        if (self.X() == v.X() and
            self.Y() == v.Y() and
            self.Z() == v.Z() and
            self.T() == v.T()):
            return True
        return False


    def __getitem__(self, key):
        if (key == 0 or key == 'x' or key == 'px'):
            return self.X()
        if (key == 1 or key == 'y' or key == 'py'):
            return self.Y()
        if (key == 2 or key == 'z' or key == 'pz'):
            return self.Z()
        if (key == 3 or key == 't' or key == 'e'):
            return self.T()
        else:
            raise AttributeError

        
    def __repr__(self):
        return '(%.4f, %.4f, %.4f, %.4f)' % (self.X(), self.Y(),
                                             self.Z(), self.T())


    def __str__(self):
        """Represents vector as a string."""
        return '(%.4f, %.4f, %.4f, %.4f)' % (self.X(), self.Y(),
                                       self.Z(), self.T())


    def __iadd__(self, v):
        """Adds another vector to self."""
        self.SetX(self.X() + v.X())
        self.SetY(self.Y() + v.Y())
        self.SetZ(self.Z() + v.Z())
        self.SetT(self.T() + v.T())


    def Generator(self):
        for i in xrange(4):
            yield self[i]
