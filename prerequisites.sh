sudo apt-get update
sudo apt-get install g++
sudo apt-get install tar
mkdir bin
#if [[ $TRAVIS_PYTHON_VERSION == 'pypy' ]]; then travis_retry pip install numpy; fi
#sudo pip install cython
#git clone https://github.com/numpy/numpy.git numpy
#cd numpy
#pip install .
#cd .. 
#sudo apt-get install $PYTHON-dev
#sudo apt-get install $PYTHON-setuptools
#sudo apt-get install $PYTHON-numpy
#sudo easy_install$PYSUF pip
#sudo pip install cython
sudo apt-get install python python-dev libatlas-base-dev gcc gfortran g++
#wget https://pypi.python.org/packages/source/s/scipy/scipy-0.17.0.tar.gz#md5=5ff2971e1ce90e762c59d2cd84837224
#tar -xvf scipy-0.17.0.tar.gz
#cd scipy-0.17.0
#python setup.py install
#cd ..
#wget https://pypi.python.org/packages/source/b/bx-python/bx-python-0.7.3.tar.gz#md5=d8c50c01c9e421bae0bbdbfa00fef6e4
#tar -xvf bx-python-0.7.3.tar.gz
#cd bx-python-0.7.3
#python setup.py install
#cd ..
#git clone https://github.com/samtools/htslib.git
#sudo make -C htslib && sudo make install -C htslib
#git clone https://github.com/samtools/samtools.git
#sudo make -C samtools && sudo make install -C samtools
#sudo cp samtools/samtools /usr/local/bin 
#sudo cp samtools/misc/wgsim /usr/local/bin
#git clone https://github.com/samtools/bcftools.git
#sudo make -C bcftools && sudo make install -C bcftools
#sudo cp bcftools/bcftools /usr/local/bin
#git clone https://github.com/lh3/bwa.git
#sudo make -C bwa
#sudo cp bwa/bwa /usr/local/bin
#wget https://github.com/broadinstitute/picard/releases/download/1.131/picard-tools-1.131.zip
#unzip picard-tools-1.131.zip
#wget http://last.cbrc.jp/last-716.zip
#unzip last-716.zip
#cd last-716
#make CXXFLAGS=-O3 -C last-716 && sudo make install -C last-716
#cd ..
#wget http://gatb-tools.gforge.inria.fr/versions/bin/minia-2.0.3-Linux.tar.gz
#tar -xvf minia-2.0.3-Linux.tar.gz
#sudo cp minia-2.0.3-Linux/bin/{dbgh5,dbginfo,h5dump,minia} /usr/local/bin
#sudo pip install cython
#sudo pip install pysam
#wget https://github.com/cython/cython/archive/master.zip
#unzip master.zip 
#cd cython-master/
#python setup.py install
#cd ..
git clone https://github.com/adamewing/align.git
cd align
python setup.py install
cd ..
