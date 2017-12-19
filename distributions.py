from __future__ import division, print_function
from random import random as _random, expovariate as _exp, gauss as _gauss
from random import uniform as _uniform, seed as _seed
from math import sqrt as _sqrt, log as _log
try:
    import pypdt
except ImportError:
    print('You will not be able to use the massPDG() function without' +
          'installing the pypdt module from pypi.')
    
_seed()


def uniform(minimum, maximum):
    """Identical to uniform found in random module."""
    return _uniform(minimum, maximum)


def randExp(lambd=250, minimum=50):
    """Return a random exponentially distributed value. Default has a 
    minimum of 50 and a stdev of 250."""
    return _exp(1/lambd) + minimum


def zMass(mZ=91.188, wZ=2.495):
    """Return a gaussian distribution. Default has a mean of 91.188 and a 
    width of 2.495."""
    return _gauss(mZ, wZ * (2*_sqrt(2*_log(2))))


def randMom(m):
    """Return a random number distributed uniformly on the interval (-m, m)."""
    return  m * (2 * _random() - 1) 


def massPDG(pdgid):
    """Gaussian mass distribution with mass and width of a chosen particle.
    Accept the Particle Data Group ID as the only argument."""
    return gauss(pypdt.get(pdgid).mass,
                 pypdt.get(pdgid).width / (2*_sqrt(2*_log(2))))

def random():
    """Return a random number on the interval [0, 1). 
    Identical to random() from the random module."""
    return _random()
    
