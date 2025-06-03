# Start from a base image that includes Miniconda
FROM continuumio/miniconda3

# Set the working directory
WORKDIR /app

# Copy the environment file into the container
COPY medaka.yaml .

# Create the Conda environment from the environment file
RUN conda env create -f medaka.yaml

# Ensure the environment is correctly configured
ENV PATH /opt/conda/envs/medaka/bin:$PATH

# Run medaka tools download_models to get all models
RUN medaka tools download_models

RUN echo "source activate medaka" > ~/.bashrc

# Set the default shell to bash
SHELL ["/bin/bash", "-c"]
