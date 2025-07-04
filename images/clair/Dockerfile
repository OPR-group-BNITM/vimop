# Start from a base image that includes Miniconda
FROM continuumio/miniconda3

# Set the working directory
WORKDIR /app

# Copy the environment file into the container
COPY clair.yaml .

# Create the Conda environment from the environment file
RUN conda env create -f clair.yaml

# Ensure the environment is correctly configured
ENV PATH /opt/conda/envs/clair/bin:$PATH

# clone rerio repository
RUN git clone https://github.com/nanoporetech/rerio.git

# download rerio models
RUN python3 rerio/download_model.py --clair3

# Set the active environment
SHELL ["/bin/bash", "-c"]
RUN echo "source activate clair" > ~/.bashrc
