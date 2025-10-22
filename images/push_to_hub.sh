#!/usr/bin/env bash

images=(canu centrifuge general ingress medaka report structural_variants)
docker_user=oprgroup

if [[ $# -gt 0 ]]
then
  images=("$@")
fi

for img in "${images[@]}"
do
    version=$(cat $img/version.txt)
    echo $img:$version
    docker push ${docker_user}/${img}:$version
    echo ""
done
