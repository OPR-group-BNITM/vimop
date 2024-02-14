# Start from a base image that includes Miniconda
FROM continuumio/miniconda3

# Set the working directory
WORKDIR /app

# Copy the Conda environment file into the image
COPY general.yaml ./general.yaml

# Create the Conda environment using the general.yaml file
RUN conda env create -f general.yaml

# Set the environment variable to activate the Conda environment
# This line is optional if you prefer to activate the environment manually
RUN echo "source activate $(head -1 general.yaml | cut -d' ' -f2)" > ~/.bashrc

# Update the PATH to include the Conda environment's binaries
ENV PATH /opt/conda/envs/$(head -1 general.yaml | cut -d' ' -f2)/bin:$PATH
