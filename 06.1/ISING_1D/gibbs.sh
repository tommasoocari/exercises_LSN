#!/bin/bash

./clean.sh

rm results*

for ((temp=5; temp<=30; temp++));
do

    cd inputs
    ./h_null.sh
    cd ..
    
    cp inputs/input$temp.dat input.dat

    sed -i '5d' input.dat
    sed -i '4 a 0' input.dat

    ./Monte_Carlo_ISING_1D.exe

    sed -i '6d' input.dat
    sed -i '5 a 1' input.dat

    ./Monte_Carlo_ISING_1D.exe
    
done

for ((temp=5; temp<=30; temp++));
do

    cd inputs
    ./h_positive.sh
    cd ..
    
    cp inputs/input$temp.dat input.dat

    sed -i '5d' input.dat
    sed -i '4 a 0' input.dat

    ./Monte_Carlo_ISING_1D.exe

    sed -i '6d' input.dat
    sed -i '5 a 1' input.dat

    ./Monte_Carlo_ISING_1D.exe
    
done

rm final.gibbs*

awk 'NR%2==1' results.ene.csv > final.gibbs.ene.csv
awk 'NR%2==1' results.mag.csv > final.gibbs.mag.csv
awk 'NR%2==1' results.chi.csv > final.gibbs.chi.csv
awk 'NR%2==1' results.heat.csv > final.gibbs.heat.csv

rm results*

