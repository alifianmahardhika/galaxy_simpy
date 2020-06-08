from vpython import vector, color, sqrt, sphere, rate, scene, gcurve,graph
from math import fsum
from random import gauss


# CONSTANTS
#====================================
# Universal gravitational constant
G = 6.673e-11 #m**3kg^(-1)s(-2)
#alpha_0
alpha = 1e-7
#Parsec unity
kparsec = 3.086e19 #m
# Solar mass in kg (assume average stellar mass)
SOLAR_MASS = 2.000e30 #kg 

#====================================
# GALACTIC PARAMETERS

MILKY_WAY_GALAXY_THICKNESS =  0.6 * kparsec

MAX_ORBITAL_RADIUS = 10 * kparsec #wikipedia
MIN_ORBITAL_RADIUS = 0.15 * kparsec

# Precalculated bounds to solar mass
MIN_SOLAR_MASS = 0.15 * SOLAR_MASS
MAX_SOLAR_MASS = 250 * SOLAR_MASS
AVG_SOLAR_MASS = 3.0 * SOLAR_MASS #source : https://arxiv.org/pdf/1306.4013v2.pdf & https://archive.briankoberlein.com/2014/10/12/galactic-scale/index.html
#==================================

#PARTICLE PARAMETRIZATION

# Milky Way contains about 300 billion stars
NUM_STARS_MILKY_WAY = 300 #wikipedia

# Graphical constants of particles
STAR_RADIUS = 0.025
#dt = 1e17
dt = 1e16 #ly in meters
#==================================

#VISUALIZATION of output

scene.width = 1300 #visualizattion
scene.height = 750
#===================================

# FUNCTIONS

# Limit x between lower and upper
def clamp(x, lower, upper):
    return max(min(x, upper), lower)

# Return the acceleration due to gravity on an object.
def g_accel(mass, radius):
    # Limit minimum radius to avoid flinging out too many particles
    radius = max(radius, MIN_ORBITAL_RADIUS)
    return G * mass / (radius **2 )

# Calculate acceleration on an object caused by galaxy
def accel(obj, galaxy):
    r_galaxy = galaxy.pos - obj.pos
    # We have a = F / m = G * m_center / r ^2
    return r_galaxy.norm() * g_accel(galaxy.mass, r_galaxy.mag)