# Start from a base image that includes Miniconda
FROM continuumio/miniconda3

# Set the working directory
WORKDIR /app

# Copy the environment.yml file into the container
COPY canu.yaml .

# Create the Conda environment from the environment.yml file
RUN conda env create -f canu.yaml

# Set the active environment
SHELL ["/bin/bash", "-c"]
RUN echo "source activate canu" > ~/.bashrc
ENV PATH /opt/conda/envs/canu/bin:$PATH
