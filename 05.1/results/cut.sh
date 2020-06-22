#!/bin/bash

for file in $(ls)
do
    echo ${file}
    sed -i '10001,$d' ${file} 
done
