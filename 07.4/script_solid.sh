#!/bin/bash

cd ../MCNVT/MonteCarlo_NVT

./clean.sh

cp config.fcc config.0
cp input.solid input.dat

sed -i '6d' input.dat
sed -i '5 a 10' input.dat
sed -i '7d' input.dat
sed -i '6 a 500' input.dat

./Monte_Carlo_NVT.exe

cp config.final config.0

./clean.sh

sed -i '6d' input.dat
sed -i '5 a 100' input.dat

./Monte_Carlo_NVT.exe

cp output.epot.0 ../../07.4/epot.solid.dat
cp output.pres.0 ../../07.4/pres.solid.dat
cp output.gave.0 ../../07.4/gave.solid.dat
