# Start from a base image with a recent stable version of Linux
FROM ubuntu:latest

# Install dependencies
RUN apt-get update && apt-get install -y wget bzip2

# Download and install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda \
    && rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH="/opt/conda/bin:$PATH"

# Copy the Conda environment file into the container
# COPY nanoflow-environment.txt .

# Create the Conda environment
RUN conda create -n nanoflow_in_docker python=3.9
RUN conda install -n nanoflow_in_docker -c anaconda pandas=1.4.2
RUN conda install -n nanoflow_in_docker -c conda-forge biopython=1.79

# Set the active environment
SHELL ["/bin/bash", "-c"]
RUN echo "source activate nanoflow_in_docker" > ~/.bashrc
ENV PATH /opt/conda/envs/nanoflow_in_docker/bin:$PATH

# Copy the Python script into the container
COPY . /app/

# Define the default command for the container
CMD python nanoflow-cli -r $run_id -i $datapath -l ${initials}_${today} --trim 2>&1 | tee -a $logfile


# # syntax=docker/dockerfile:1

# # Comments are provided throughout this file to help you get started.
# # If you need more help, visit the Dockerfile reference guide at
# # https://docs.docker.com/go/dockerfile-reference/

# # Want to help us make this template better? Share your feedback here: https://forms.gle/ybq9Krt8jtBL3iCk7

# ################################################################################
# # Pick a base image to serve as the foundation for the other build stages in
# # this file.
# #
# # For illustrative purposes, the following FROM command
# # is using the alpine image (see https://hub.docker.com/_/alpine).
# # By specifying the "latest" tag, it will also use whatever happens to be the
# # most recent version of that image when you build your Dockerfile.
# # If reproducability is important, consider using a versioned tag
# # (e.g., alpine:3.17.2) or SHA (e.g., alpine@sha256:c41ab5c992deb4fe7e5da09f67a8804a46bd0592bfdf0b1847dde0e0889d2bff).
# FROM alpine:latest as base

# ################################################################################
# # Create a stage for building/compiling the application.
# #
# # The following commands will leverage the "base" stage above to generate
# # a "hello world" script and make it executable, but for a real application, you
# # would issue a RUN command for your application's build process to generate the
# # executable. For language-specific examples, take a look at the Dockerfiles in
# # the Awesome Compose repository: https://github.com/docker/awesome-compose
# FROM base as build
# RUN echo -e '#!/bin/sh\n\
# echo Hello world from $(whoami)! In order to get your application running in a container, take a look at the comments in the Dockerfile to get started.'\
# > /bin/hello.sh
# RUN chmod +x /bin/hello.sh

# ################################################################################
# # Create a final stage for running your application.
# #
# # The following commands copy the output from the "build" stage above and tell
# # the container runtime to execute it when the image is run. Ideally this stage
# # contains the minimal runtime dependencies for the application as to produce
# # the smallest image possible. This often means using a different and smaller
# # image than the one used for building the application, but for illustrative
# # purposes the "base" image is used here.
# FROM base AS final

# # Create a non-privileged user that the app will run under.
# # See https://docs.docker.com/go/dockerfile-user-best-practices/
# ARG UID=10001
# RUN adduser \
#     --disabled-password \
#     --gecos "" \
#     --home "/nonexistent" \
#     --shell "/sbin/nologin" \
#     --no-create-home \
#     --uid "${UID}" \
#     appuser
# USER appuser

# # Copy the executable from the "build" stage.
# COPY --from=build /bin/hello.sh /bin/

# # What the container should run when it is started.
# ENTRYPOINT [ "/bin/hello.sh" ]
