#!/bin/bash

./clean.sh

cp coord_circ.dat input.coordinates

./TSP.exe

cp output.ave.0 results.circ.ave
cp output.min.0 results.circ.min
cp output.path.0 results.circ.path




