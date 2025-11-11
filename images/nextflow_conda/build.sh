#!/usr/bin/env bash

version=$(cat version.txt)
echo $version
docker build --no-cache -t oprgroup/nextflow_conda:$version ./
