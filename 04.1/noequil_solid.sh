#!/bin/bash

PHASE=solid

cd ../MolecularDynamics_NVE

./clean.sh

cp config.fcc config.0
cp input.${PHASE} input.dat

./MolDyn_NVE.exe

cp all_epot.dat ../04.1/equil_${PHASE}_epot.dat
cp all_etot.dat ../04.1/equil_${PHASE}_etot.dat
cp all_ekin.dat ../04.1/equil_${PHASE}_ekin.dat
cp all_temp.dat ../04.1/equil_${PHASE}_temp.dat
