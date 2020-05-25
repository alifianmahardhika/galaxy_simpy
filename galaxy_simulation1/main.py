from globals import *
from galaxy import Galaxy
import sys

def main():

    t = 0
    milky_way = Galaxy(
        num_stars=NUM_STARS_MILKY_WAY,
        pos=vector(-5, 0, 0) * kparsec,
        vel=vector(0, 0, 0),
        radius=MAX_ORBITAL_RADIUS,
        thickness=MILKY_WAY_GALAXY_THICKNESS,
        color=vector(0.9, 0.9, 1)
    )
    posgraph = gcurve(color=color.green)
    

    while True:
        rate(100) #part of visualization in vpython
        

        for i in range(len(milky_way.stars)):
            star = milky_way.stars[i]
            star.pos += star.vel * dt
            star.vel += accel(star, milky_way) * dt
                
        t += dt
        
if __name__ == '__main__':
    main()
    
