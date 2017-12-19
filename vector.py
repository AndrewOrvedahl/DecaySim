from math import sqrt as _sqrt, atan2 as _atan2, log as _log, pi as _pi


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
        """Return the squared momentum."""
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
        return (self._x * v._x + self._y * v._y + self._z * v._z)


    def Unit3(self):
        """Return a unit three-vector."""
        mag = self.P()
        return VecThree(self._x / mag,
                        self._y / mag,
                        self._z / mag)

    
    
class VecFour(VecThree):
    """Simple four-vector implementation. Every method for retrieving or
    setting a value is implemented for both momentum and position. The 
    methods are interchangeable."""
    def __init__(self, x=0, y=0, z=0, t=0):
        super(VecFour, self).__init__(x, y, z)
        self._t = t
        return


    def M2(self):
        """Return E^2 - P^2."""
        return self._t**2 - (self._x**2 + self._y**2 + self._z**2)


    def T(self):
        """Return the t-component."""
        return self._t


    def E(self):
        """Return the e-component."""
        return self._t
    

    def SetT(self, t):
        """Set the t-component."""
        self._t = t


    def SetE(self, e):
        """Set the e-component."""
        self._t = e


    def Dot4(self, v):
        """Return E1*E2 - P1*P2."""
        return self._t * v._t - self.Dot3(v)


    def Unit4(self):
        """Return a unit four-vector."""
        mag = _sqrt(self._x**2 + self._y**2 + self._z**2 + self._t**2)
        return VecFour(self._x / mag,
                       self._y / mag,
                       self._z / mag,
                       self._t / mag)

    
