# Start from a base image that includes Miniconda
FROM continuumio/miniconda3

# Set the working directory
WORKDIR /app

# Copy the environment.yml file into the container
COPY centrifuge.yaml .

# Create the Conda environment from the environment.yml file
RUN conda env create -f centrifuge.yaml

# Set the active environment
SHELL ["/bin/bash", "-c"]
RUN echo "source activate centrifuge" > ~/.bashrc
ENV PATH /opt/conda/envs/centrifuge/bin:$PATH

