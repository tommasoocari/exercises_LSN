#!/bin/bash

PHASE=gas

./equil_${PHASE}.sh

./clean.sh

./MolDyn_NVE.exe

cp output.epot.0 results/${PHASE}_epot.out
cp output.etot.0 results/${PHASE}_etot.out
cp output.ekin.0 results/${PHASE}_ekin.out
cp output.temp.0 results/${PHASE}_temp.out
cp output.gave.0 results/${PHASE}_gave.out
cp output.gofr.0 results/${PHASE}_gofr.out
