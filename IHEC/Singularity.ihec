BootStrap: docker
From: ubuntu:xenial

%runscript
    exec /usr/local/bin/gemBS $@

%help
    gemBS singularity container
	 
%post
	 (mkdir /ext && cd /ext && mkdir disk1 disk2 disk3 disk4 disk5 disk6 disk7 disk8 disk9)
    apt-get update
	 apt-get install -y libpng-dev uuid-dev libmysqlclient-dev
	 apt-get install -y python3 build-essential git python3-pip wget pigz
    apt-get install -y zlib1g-dev libbz2-dev gsl-bin libgsl0-dev
    apt-get install -y libncurses5-dev liblzma-dev libssl-dev libcurl4-openssl-dev
    pip3 install 'matplotlib<3.0' multiprocess
    mkdir /usr/local/build; cd /usr/local/build
	 git clone --recursive https://github.com/heathsc/gemBS.git
    (cd gemBS; python3 setup.py install)
    rm -rf gemBS && cd && rmdir /usr/local/build
