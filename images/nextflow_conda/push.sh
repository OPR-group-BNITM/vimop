#!/usr/bin/env bash

version=$(cat version.txt)
echo $version
docker push oprgroup/nextflow_conda:$version
