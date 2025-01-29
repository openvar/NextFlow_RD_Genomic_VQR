# Use the continuumio/miniconda3 base image
FROM continuumio/miniconda3:23.3.1-0

# Update the package list and install necessary packages using apt
# Note: I'm running a Mac M1, and these packages do not seem to be available
# for ARM via conda, so using apt to install
RUN conda init bash && \
    echo "conda activate base" > ~/.bashrc && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda install -y gatk4==4.3.0.0 && \
    conda clean --all -f -y && \
    apt-get update && \
    apt-get install -y samtools=1.11-1 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    apt-get remove -y && \
    apt-get autoremove -y

# Metadata indicating the author or maintainer of the image
LABEL maintainer="Peter J Freeman <peter.j.freeman@manchester.ac.uk>"

# Metadata about the purpose or description of the image
LABEL description="Docker image containing gatk4==4.3.0.0"

# Default command to run when the container starts
CMD ["/bin/bash", "-l"]