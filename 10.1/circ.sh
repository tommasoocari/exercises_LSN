#!/bin/bash

./clean.sh

cp coord_circ.dat input.coordinates

./TSP.exe

cp output.fitness.0 results.circ.fitness
cp output.path.0 results.circ.path




