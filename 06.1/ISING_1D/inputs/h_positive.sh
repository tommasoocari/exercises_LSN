#!/bin/bash

for ((temp=5; temp<=30; temp++));
do
    sed -i '4d' input$temp.dat
    sed -i '3 a 0.02' input$temp.dat  

done
