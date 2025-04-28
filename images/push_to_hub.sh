#!/usr/bin/env bash

images=(canu centrifuge general ingress medaka report)
images=(medaka)
docker_user=oprgroup

for img in "${images[@]}";
do
    version=$(cat $img/version.txt)
    echo $img:$version
    docker push ${docker_user}/${img}:$version
    echo ""
done
