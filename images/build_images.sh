#!/usr/bin/env bash

images=(canu centrifuge general ingress medaka report)

for img in "${images[@]}";
do
    echo $img
    docker build --no-cache -t oprgroup/${img}:0.0.1 $img/
    echo ""
done
