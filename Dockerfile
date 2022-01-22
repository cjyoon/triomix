# Triomix v1.3 2022.01.21 
FROM ubuntu:focal
ENV DEBIAN_FRONTEND noninteractive
# Upgrade installed packages
RUN apt-get update && apt-get install -y -f autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev libncurses5-dev bwa bedtools build-essential python3.8 python3.8-venv python3-pip python3.8-dev git wget vim less software-properties-common python-is-python3 openjdk-8-jdk locales locales-all curl


ENV R_HOME=/usr/local/lib/R
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
    add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' && \
    apt update && apt install -y r-base-dev libxml2 libxml2-dev
# Install R packages
RUN R -e "install.packages('BiocManager', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install('DNAcopy')"
RUN R -e "BiocManager::install('aroma.light')"
RUN R -e "install.packages('PSCBS', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('tidyverse', dependencies=TRUE)"
RUN R -e "install.packages('BiocManager', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('optparse', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('bbmle', dependencies=TRUE, repos='http://cran.rstudio.com/')"

# Install Htslib
RUN mkdir tools && cd tools && wget https://github.com/samtools/htslib/releases/download/1.11/htslib-1.11.tar.bz2 && \
    tar xvfj htslib-1.11.tar.bz2 && \
    cd htslib-1.11 &&  ./configure && make && make install && cd $HOME 

# Install samtools
RUN cd tools && wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
    tar xvfj samtools-1.9.tar.bz2 && \
    cd samtools-1.9 &&  ./configure && make && make install && cp samtools /usr/local/bin/ && cd $HOME

# Install python3 packages
RUN pip3 install pysam snakemake numpy pandas scipy 

ENV PATH="/usr/local/bin/:${PATH}"

# Clean up
RUN cd tools && rm -f *.tar.bz2 && rm -f *.tar.gz && cd $HOME 

RUN cd tools && git clone https://github.com/cjyoon/triomix.git 
