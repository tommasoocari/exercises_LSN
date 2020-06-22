#!/bin/bash

./clean.sh

cp coord_square.dat input.coordinates

mpiexec -np 4 TSP.exe

for i in 0 1 2 3
do
    cp output.ave.$i results/results.square.ave.$i
    cp output.min.$i results/results.square.min.$i
    cp output.path.$i results/results.square.path.$i
done
