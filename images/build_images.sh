#!/usr/bin/env bash

set -euo pipefail

images=(canu centrifuge general ingress medaka report structural_variants)

if [[ $# -gt 0 ]]
then
  images=("$@")
fi

for img in "${images[@]}";
do
    version=$(cat $img/version.txt)
    echo $img:$version
    docker build --no-cache -t oprgroup/${img}:$version $img/
    echo ""
done
