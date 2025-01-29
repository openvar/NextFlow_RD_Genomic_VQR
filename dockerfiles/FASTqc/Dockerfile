# Use the continuumio/miniconda3 base image
FROM continuumio/miniconda3:23.3.1-0

# Update the package list and install necessary packages using apt
# Note: I'm running a Mac M1, and these packages do not seem to be available
# for ARM via conda, so using apt to install
RUN conda init bash && \
    echo "conda activate base" > ~/.bashrc && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda install -y \
        fastqc==0.12.1 && \
    conda clean --all -f -y

# Metadata indicating the author or maintainer of the image
LABEL maintainer="Peter J Freeman <peter.j.freeman@manchester.ac.uk>"

# Metadata about the purpose or description of the image
LABEL description="Docker image containing FASTqc v0.12.1"

# Default command to run when the container starts
CMD ["/bin/bash", "-l"]