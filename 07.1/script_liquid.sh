#!/bin/bash

cd ../MCNVT/MonteCarlo_NVT

./clean.sh

cp config.fcc config.0
cp input.liquid input.dat

sed -i '6d' input.dat
sed -i '5 a 3' input.dat
sed -i '7d' input.dat
sed -i '6 a 5000' input.dat

./Monte_Carlo_NVT.exe

cp config.final config.0

./clean.sh

sed -i '6d' input.dat
sed -i '5 a 100' input.dat

./Monte_Carlo_NVT.exe

cp all_epot.dat ../../07.1/all_epot_liquid.dat
cp all_pres.dat ../../07.1/all_pres_liquid.dat

cp output.gave.0 ../../07.1/gave_metropolis_liquid.dat

