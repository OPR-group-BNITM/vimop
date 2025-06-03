#!/usr/bin/env bash

images=(report)

for img in "${images[@]}";
do
    version=$(cat $img/version.txt)
    echo $img:$version
    docker build --no-cache -t oprgroup/${img}:$version $img/
    echo ""
done
