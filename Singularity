BootStrap: docker
From: ubuntu:bionic   # This is a comment

%runscript
    echo "This is what happens when you run the container..."

%post
    apt-get update
    apt-get install -y python build-essential git python-pip wget pigz
    apt-get install -y zlib1g-dev libbz2-dev gsl-bin libgsl0-dev
    apt-get install -y libncurses5-dev liblzma-dev libssl-dev libcurl4-openssl-dev
    (cd /usr/local/bin; \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig \
    && chmod 755 wigToBigWig)
    pip install numpy
    pip install matplotlib
    mkdir /usr/local/build; cd /usr/local/build
    wget https://github.com/samtools/samtools/releases/download/1.8/samtools-1.8.tar.bz2
    tar -jxf samtools-1.8.tar.bz2 && rm samtools-1.8.tar.bz2
    (cd samtools-1.8; ./configure --prefix=/usr/local && \
    make all all-htslib && make install install-htslib)
    rm -rf samtools-1.8
    wget https://github.com/samtools/bcftools/releases/download/1.8/bcftools-1.8.tar.bz2
    tar -jxf bcftools-1.8.tar.bz2 && rm bcftools-1.8.tar.bz2
    (cd bcftools-1.8; ./configure --prefix=/usr/local && make && make install)
    rm -rf bcftools-1.8
	 git clone --recursive https://github.com/heathsc/gemBS.git
    (cd gemBS; python setup.py install && \
	 cp extras/gemBS_singularity.py /usr/local/bin/gemBS && chmod 755 /usr/local/bin/gemBS)
    rm -rf gemBS && cd && rmdir /usr/local/build
