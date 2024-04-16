#!/usr/bin/env bash

images=(canu centrifuge general medaka)

for img in "${images[@]}";
do
    echo $img
    docker build -t opr_${img}:0.0.1 $img/
    echo ""
done