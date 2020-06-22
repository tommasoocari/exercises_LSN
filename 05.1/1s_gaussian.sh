#!/bin/bash

cp 1s_gauss.dat input.dat
cp 1s_distribution.cpp distribution.cpp

make

./main.exe
