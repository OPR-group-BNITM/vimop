# This is not an image used by ViMOP. Instead, it is an image
# that has conda and nextflow installed. Building a container from this
# image and mounting a directory with ViMOP and the ViMOP_DB, one can run
# ViMOP only in this container. This can be used, to run ViMOP on a computing cluster,
# e.g. in a slurm system.

# ---- Base OS ----
FROM continuumio/miniconda3

# Set environment variables to avoid prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Set the working directory
WORKDIR /app

# ---- System deps (incl. Java for Nextflow) ----
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl wget ca-certificates bash bzip2 git openssh-client \
    openjdk-17-jre-headless \
    && rm -rf /var/lib/apt/lists/*

# ---- Nextflow ----
# Optionally pin a Nextflow version with: --build-arg NXF_VER=24.10.3
ARG NXF_VER
RUN if [ -n "$NXF_VER" ]; then export NXF_VER="$NXF_VER"; fi \
    && curl -fsSL https://get.nextflow.io | bash \
    && mv nextflow /usr/local/bin/ \
    && chmod +x /usr/local/bin/nextflow
