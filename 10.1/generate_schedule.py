import numpy as np

file_out = open("schedule.dat","w")

for i in range(1000):
    file_out.write("1000 " + str(i) + "\n")

file_out.close()
    
