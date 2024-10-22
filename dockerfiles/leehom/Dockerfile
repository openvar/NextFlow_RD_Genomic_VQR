# Use the official Ubuntu base image
FROM ubuntu:20.04

# Install necessary packages and clone Leehom from GitHub
RUN apt-get update && \
    apt-get install -y \
    git \
    cmake \
    g++ \
    make \
    zlib1g-dev && \
    git clone https://github.com/nextflow-io/leehom.git /leehom && \
    cd /leehom && \
    make && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Metadata indicating the author or maintainer of the image
LABEL maintainer="Peter J Freeman <peter.j.freeman@manchester.ac.uk>"

# Metadata about the purpose or description of the image
LABEL description="Docker image for building Leehom for cross-platform usage"

# Set the default command to run when the container starts
CMD ["/bin/bash", "-l"]
