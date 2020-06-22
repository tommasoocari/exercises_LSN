#!/bin/bash

./clean.sh

cp coord_circ.dat input.coordinates

mpiexec -np 4 TSP.exe

for i in 0 1 2 3
do
    cp output.ave.$i results/results.circ.ave.$i
    cp output.min.$i results/results.circ.min.$i
    cp output.path.$i results/results.circ.path.$i
done



