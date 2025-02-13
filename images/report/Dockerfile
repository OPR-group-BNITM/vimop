# Start from a base image that includes Miniconda
FROM continuumio/miniconda3

# Set the working directory
WORKDIR /app

# Copy your repository into the container
COPY . .

# Create the Conda environment from the environment.yml file
RUN conda env create -f report.yaml && \
    echo "source activate report" > ~/.bashrc


# Navigate to the `ezcharts` directory and install the requirements
RUN /bin/bash -c "source activate report && cd ezcharts_opr && pip install -r requirements.txt && python setup.py install"

# Set the active environment
SHELL ["/bin/bash", "-c"]
ENV PATH /opt/conda/envs/report/bin:$PATH
