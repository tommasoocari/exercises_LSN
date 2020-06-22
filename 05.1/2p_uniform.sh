#!/bin/bash

cp 2p_uniform.dat input.dat
cp 2p_distribution.cpp distribution.cpp

make

./main.exe
