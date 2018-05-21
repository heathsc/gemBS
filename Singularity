BootStrap: docker
From: ubuntu:bionic   # This is a comment

%runscript
    exec /usr/local/bin/gemBS $@

%help
    gemBS singularity container
	 
%post
    apt-get update
    apt-get install -y python build-essential git python-pip wget pigz
    apt-get install -y zlib1g-dev libbz2-dev gsl-bin libgsl0-dev
    apt-get install -y libncurses5-dev liblzma-dev libssl-dev libcurl4-openssl-dev
    (cd /usr/local/bin; \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig \
    && chmod 755 wigToBigWig)
    pip install numpy matplotlib configparser multiprocess
    mkdir /usr/local/build; cd /usr/local/build
	 (mkdir /ext && cd /ext && mkdir disk1 disk2 disk3 disk4 disk5 disk6 disk7 disk8 disk9)
	 git clone --recursive https://github.com/heathsc/gemBS.git -b development
    (cd gemBS; python setup.py install)
    rm -rf gemBS && cd && rmdir /usr/local/build
