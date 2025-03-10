#!/usr/bin/env bash

images=(canu centrifuge general ingress medaka report)
docker_user=oprgroup

for img in "${images[@]}";
do
    echo $img
    docker push ${docker_user}/${img}:0.0.1
    echo ""
done
