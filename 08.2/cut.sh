#!/bin/bash

for file in $(ls)
do
    echo ${file}
    sed -i '1000001,$d' ${file} 
done
