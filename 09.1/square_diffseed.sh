#!/bin/bash

./clean.sh

cp seed2.in seed.in

cp coord_square.dat input.coordinates

./TSP.exe

cp output.ave.0 results.square.ave_diff
cp output.min.0 results.square.min_diff
cp output.path.0 results.square.path_diff

cp seed1.in seed.in
