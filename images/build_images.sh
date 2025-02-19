#!/usr/bin/env bash

images=(canu centrifuge general ingress medaka report)

for img in "${images[@]}";
do
    echo $img
    docker build --no-cache -t opr_${img}:0.0.1 $img/
    echo ""
done
