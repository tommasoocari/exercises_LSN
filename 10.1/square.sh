#!/bin/bash

./clean.sh

cp coord_square.dat input.coordinates

./TSP.exe

cp output.fitness.0 results.square.fitness
cp output.path.0 results.square.path
