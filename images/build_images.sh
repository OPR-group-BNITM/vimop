#!/usr/bin/env bash

images=(canu centrifuge general ingress medaka report structural_variants clair)

for img in "${images[@]}";
do
    version=$(cat $img/version.txt)
    echo $img:$version
    docker build --no-cache -t oprgroup/${img}:$version $img/
    echo ""
done
