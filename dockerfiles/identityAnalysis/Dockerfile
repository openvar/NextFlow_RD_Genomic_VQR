# Use the continuumio/miniconda3 base image
FROM continuumio/miniconda3:23.3.1-0

# Update the package list and install necessary packages using conda
RUN conda init bash && \
    echo "conda activate base" > ~/.bashrc && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda install -y haplogrep==2.4.0 && \
    conda install -y bcftools && \
    apt-get update && \
    apt-get install -y unzip && \
    apt-get clean && \
    apt-get install -y build-essential && \
    apt-get install -y libz-dev && \
    apt-get install -y wget && \
    apt-get install -y cmake && \
    apt-get install -y liblzma-dev && \
    apt-get install -y libbz2-dev && \
    apt-get install -y libcurl4-openssl-dev && \
    apt-get install -y zlib1g-dev && \
    rm -rf /var/lib/apt/lists/* && \
    apt-get remove -y && \
    apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/* && \
    conda clean --all -f -y

# Install HTSlib
RUN wget https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2 && \
    tar -xvjf htslib-1.21.tar.bz2 && \
    cd htslib-1.21 && \
    make && \
    make install && \
    cd .. && \
    rm -rf htslib-1.21 htslib-1.21.tar.bz2

# Download and unzip plink2
RUN wget https://s3.amazonaws.com/plink2-assets/alpha5/plink2_linux_x86_64_20241020.zip && \
    unzip plink2_linux_x86_64_20241020.zip && \
    mv plink2 /usr/local/bin/ && \
    chmod +x /usr/local/bin/plink2 && \
    rm plink2_linux_x86_64_20241020.zip

# Metadata indicating the author or maintainer of the image
LABEL maintainer="Peter J Freeman <peter.j.freeman@manchester.ac.uk>"

# Metadata about the purpose or description of the image
LABEL description="Docker image containing bcftools haplogrep2 and plink for DNA identity analysis."

# Default command to run when the container starts
CMD ["/bin/bash", "-l"]
