from vpython import vector, color, sqrt, sphere, rate, scene, gcurve,graph
from math import fsum
from random import gauss
import numpy as np

# CONSTANTS

# Universal gravitational constant
G = 6.673e-11 #m**3kg^(-1)s(-2)

#alpha_0
alpha = 1e-7

scene.width = 1300 #visualizattion
scene.height = 750

# Solar mass in kg (assume average stellar mass)
SOLAR_MASS = 2.000e30 #kg 

# Precalculated bounds to solar mass
MIN_SOLAR_MASS = SOLAR_MASS * 0.15
MAX_SOLAR_MASS = SOLAR_MASS * 250
AVG_SOLAR_MASS = SOLAR_MASS * 2.0

#Parsec unity
parsec = 3.086e19 #m

# Galactic parameters
MAX_ORBITAL_RADIUS = parsec * 10 #wikipedia
MIN_ORBITAL_RADIUS = parsec * 0.15

MILKY_WAY_GALAXY_THICKNESS = parsec * 0.6 #wikipedia
#ANDROMEDA_GALAXY_THICKNESS = DIST_SCALE * 0.2


# Milky Way contains about 300 billion stars
NUM_STARS_MILKY_WAY = 700

# Graphical constants of particles
STAR_RADIUS = 0.025
#dt = 1e17
dt = 1e16

# FUNCTIONS

# Limit x between lower and upper
def clamp(x, lower, upper):
    return max(min(x, upper), lower)


# Return the force due to gravity on an object (Newton universal gravity law)
def gravity(mass1, mass2, radius):
    return G * mass1 * mass2 / radius**2

# Return the acceleration due to gravity on an object.
def g_accel(mass, radius):
    # Limit minimum radius to avoid flinging out too many particles
    radius = max(radius, MIN_ORBITAL_RADIUS)
    return G * mass / radius / radius

# Calculate acceleration on an object caused by galaxy
def accel(obj, galaxy):
    r_galaxy = galaxy.pos - obj.pos
    # We have a = F / m = G * m_center / r ^2
    return r_galaxy.norm() * g_accel(galaxy.mass, r_galaxy.mag)
