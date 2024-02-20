# Start from a base image that includes Miniconda
FROM continuumio/miniconda3

# Set the working directory
WORKDIR /app

# Copy the environment.yml file into the container
COPY general.yaml .

# Create the Conda environment from the environment.yml file
RUN conda env create -f general.yaml

# Set the active environment
SHELL ["/bin/bash", "-c"]
RUN echo "source activate general" > ~/.bashrc
ENV PATH /opt/conda/envs/general/bin:$PATH

# Set the default command to run when the container starts
CMD ["samtools", "--version"]
