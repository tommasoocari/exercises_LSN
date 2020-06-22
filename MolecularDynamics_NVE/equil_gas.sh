#!/bin/bash

PHASE=gas

./clean.sh

cp config.fcc config.0
cp input.${PHASE} input.dat

sed -i '1d' input.dat
sed -i '1i1.0' input.dat

sed -i '9d' input.dat
sed -i '8 a 10' input.dat 

./MolDyn_NVE.exe

cp old.final old.0
cp config.final config.0

sed -i '6d' input.dat
sed -i '5 a 1' input.dat

sed -i '7d' input.dat
sed -i '6 a 1' input.dat

sed -i '8d' input.dat
sed -i '7 a 1.175' input.dat 

sed -i '9d' input.dat
sed -i '8 a 100' input.dat 

