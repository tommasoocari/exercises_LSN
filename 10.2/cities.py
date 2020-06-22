import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math
import random as rnd

rnd.seed(0)

circ = []

def uniform_angle():
    x = rnd.uniform(-1,1)
    y = rnd.uniform(-1,1)
    r = math.sqrt(x*x + y*y)
    while(x*x + y*y < 1.0):
        x = rnd.uniform(-1,1)
        y = rnd.uniform(-1,1)
        r = math.sqrt(x*x + y*y)
    return [x / r,y / r]    

for i in range(32):
    circ.append(uniform_angle())

circ = np.array(circ)

f = open("coord_circ.dat","w")

for i in circ:
    f.write(str(i[0])+ " " +str(i[1]) + "\n")
    
    
f.close()

square = []

for i in range(32):
    square.append([rnd.uniform(0,1), rnd.uniform(0,1)])

square = np.array(square)

f = open("coord_square.dat","w")

for i in square:
    f.write(str(i[0])+ " " +str(i[1]) + "\n")
    
f.close()
