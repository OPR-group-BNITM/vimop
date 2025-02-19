# Use an official Ubuntu as a base image
FROM ubuntu:20.04

# Set environment variables to avoid prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies
RUN apt-get update && \
    apt-get install -y \
    build-essential \
    g++ \
    wget \
    zlib1g-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    openjdk-8-jdk \
    git \
    && rm -rf /var/lib/apt/lists/*

# Set Canu version
ENV CANU_VERSION=2.2

# Download and install Canu
RUN wget https://github.com/marbl/canu/releases/download/v$CANU_VERSION/canu-$CANU_VERSION.Linux-amd64.tar.xz && \
    tar -xvf canu-$CANU_VERSION.Linux-amd64.tar.xz && \
    mv canu-$CANU_VERSION /usr/local/canu && \
    chmod +x /usr/local/canu/bin/canu && \
    ln -s /usr/local/canu/bin/canu /usr/local/bin/canu

# Add Canu to PATH
ENV PATH="/usr/local/canu/bin:$PATH"