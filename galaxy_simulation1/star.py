from globals import *

class Star(object):
    def __init__(self, mass, radius, pos, vel, color):
        self.obj = sphere(pos=pos / parsec, radius=radius, color=color)
        self.mass = mass
        self.vel = vel
        self._pos = pos

    # Externally use scaled version for physics, use normalized version for graphics
    @property
    def pos(self):
        return self._pos

    @pos.setter
    def pos(self, value):
        self.obj.pos = value / parsec
        self._pos = value

    def __str__(self):
        return "Mass: " + str(self.mass) + "\nPos: " + str(self.pos) + \
            "\nVel: " + str(self.vel)
