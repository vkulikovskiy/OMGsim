OMGsim documentation {#mainpage}
----------------------------------
by Christophe Hugon - CPPM (CNRS) and Vladimir Kulikovskiy - INFN Genova

##Introduction

OMGsim software objective is to provide an easy to use simulation of the neutrino telescopes in general, 
ANTARES, NEMO and KM3NeT in particular. An important attention as been given to the simulation of the 
photocathode, using a dedicated thin layer and complex index simulation.

##Features
This simulation will allow you to:

- easily determine a precise geometry of the detector
- easily determine a pricise geometry of the phototubes and optical modules
- easily determine the different materials and properties of the detector and phototubes
- use the exact Petzold water scattering of light (with the base parameters measured for ANTARES)
- Test any experiment you would like to reproduce
     - Laboratory angular acceptance
     - 40K decay in water
     -...
     
##Description
- Quick start [here](md_documentation_pages_quickstart.html).
- KM3NeT DOM Simulation [here](md_documentation_pages_DOM.html).
- Calibration [here](md_documentation_pages_calibration.html).
- Visualization [here] (md_documentation_pages_visualization.html).
- Output file format [here] (md_documentation_pages_outputformat.html).


##Installation

###Warning
This tool uses some c++11 great features as auto declaration and in-header attribute initialization. 
It is strongly advised to have a newer than 4.6 gcc version (even more recent advised).

###Dependencies
The dependencies are already satisfied at the CC de Lyon. 
So, except if you strongly want to install it on your own computer, goes directly to the how to at Lyon.

- GCC >= 4.7
- cmake >= 3.2
- ROOT > 5.32
- libconfig
- GEANT4 > 10
- Doxygen > 1.8

###Download and compilation
###Lyon configuration
#### SL6
The environment was prepared at Lyon to allow an easy installation. You need to set the following variables in your .bashrc or .zshrc:

    #Geant4
    cd /sps/km3net/users/hugon/sw/geant4.10.01.p01/bin
    source geant4.sh
    cd -

    #Libconfig
    export LIBCONFIG_ROOT=/sps/km3net/users/hugon/sw/libconfig-1.5
    export LD_LIBRARY_PATH=${LIBCONFIG_ROOT}/lib:${LD_LIBRARY_PATH}
    export CPLUS_INCLUDE_PATH=${LIBCONFIG_ROOT}/include/

    #GCC 5.10
    export LD_LIBRARY_PATH=/sps/km3net/users/hugon/sw/gcc/gcc-5.1.0/lib:/sps/km3net/users/hugon/sw/gcc/gcc-5.1.0/lib64:${LD_LIBRARY_PATH}
    export PATH=/sps/km3net/users/hugon/sw/gcc/gcc-5.1.0/bin:${PATH}

    #ROOT
    export ROOTSYS=/usr/local/root/v5.34.23
    cd ${ROOTSYS}
    source bin/thisroot.sh
    cd -

This OMGsim can be obtained via svn, using the command

    git clone https://github.com/vkulikovskiy/OMGsim.git
    cd trunk

then to compile it

    mkdir build
    cd build
    cmake .. -DCMAKE_CXX_COMPILER=g++-510 -DCMAKE_C_COMPILER=gcc-510
    make


#### CentOS7
A common installation of Geant4 and some libraries was prepared, so the configuration should be the following:

    #Geant4
    cd $KM3NET_THRONG_DIR/share/geant/v10.01.p03/bin/
    source geant4.sh
    cd -

    #Libconfig
    export LIBCONFIG_ROOT=/sps/km3net/users/hugon/sw/libconfig-1.5
    export LD_LIBRARY_PATH=${LIBCONFIG_ROOT}/lib:${LD_LIBRARY_PATH}
    export CPLUS_INCLUDE_PATH=${LIBCONFIG_ROOT}/include

    #GCC 5.10
    export GCC_ROOT=/pbs/software/x86_64_el7/gcc/5.5.0
    export LD_LIBRARY_PATH=${GCC_ROOT}/lib:${GCC_ROOT}/lib64:${LD_LIBRARY_PATH}
    #export PATH=${GCC_ROOT}/bin:/sps/km3net/users/kulikovs/OMGsim/gcc:${PATH} 
    export PATH=${GCC_ROOT}/bin:${PATH}

    #ROOT
    export ROOTSYS=/usr/local/root/v5.34.23
    cd ${ROOTSYS}
    source bin/thisroot.sh
    cd -

This OMGsim can be obtained via svn, using the command

    git clone https://github.com/vkulikovskiy/OMGsim.git
    cd trunk

Then to compile it:

    #To force usage GCC 5.5.0 compiler (note that this trick was also used to compile geant4 in $KM3NET_THRONG_DIR/share)
    export CC=${GCC_ROOT}/bin/gcc     
    export CXX=${GCC_ROOT}/bin/g++
    mkdir build
    cd build
    cmake ..
    make


### MacOS installation 
Installation of missing packages is faster with HomeBrew packet manager.
Install it with:

    /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

Install Motif libraries:

    brew install motif

If cmake is missing, one could install it with HomeBrew as well:

    brew install cmake

Download and compile geant4.10.01.p03:
    
    wget http://geant4.web.cern.ch/geant4/support/source/geant4.10.01.p03.tar.gz
    tar xf geant4.10.01.p03.tar.gz
    mkdir geant4.10.01.p03-build
    cd geant4.10.01.p03-build
    cmake -DCMAKE_INSTALL_PREFIX=../geant4.10.01.p03-install ../geant4.10.01.p03 -DGEANT4_INSTALL_DATA=ON -DGEANT4_USE_XM=ON
    make -j4 install

Add to your ~/.bashrc  or ~/.zshrc

    cd ~/Soft/geant4.10.01.p03-install/bin/
    source geant4.sh
    cd -

#### Root installation
Download and compile ROOT:

    wget https://root.cern.ch/download/root_v6.10.08.source.tar.gz
    tar xf root_v6.10.08.source.tar.gz
    mkdir root_v6.10.08_build
    cd root_v6.10.08_build
    cmake ..
    make -j4

Add to your ~/.bashrc

    cd ~/Soft/root-6.10.08_build
    souce bin/thisroot.sh
    cd -

#### OMGsim
OMGsim needs libconfig++

    brew install libconfig

Install SVN:

    brew install subversion

On Mac X11 libraries are installed not in common location. Add to your ~/.bashrc or ~/.zshrc (and reload
the environment): 

    export X11HEADERS_PATH=/opt/X11/include

Download the code and compile:

    git clone https://github.com/vkulikovskiy/OMGsim.git 
    mkdir build
    cd build
    cmake ..
    make -j4

### Ubuntu installation

#### Geant4 Installation
Some libraries are needed:

    sudo apt install cmake libexpat1-dev libxerces-c-dev qt4-dev-tools libxmu-dev libmotif-dev

Download and compile geant4.10.01.p03:
    
    wget http://geant4.web.cern.ch/geant4/support/source/geant4.10.01.p03.tar.gz
    tar xf geant4.10.01.p03.tar.gz
    mkdir geant4.10.01.p03-build
    cd geant4.10.01.p03-build
    cmake -DCMAKE_INSTALL_PREFIX=../geant4.10.01.p03-install ../geant4.10.01.p03 -DGEANT4_INSTALL_DATA=ON -DGEANT4_USE_XM=ON
    make -j4 install

To configure in graphical mode

    sudo apt install cmake-curses-gui
    ccmake  -DCMAKE_INSTALL_PREFIX=../geant4.10.01.p03-install ../geant4.10.01.p03 -DGEANT4_INSTALL_DATA=ON -DGEANT4_USE_XM=ON

Add to your ~/.bashrc:

    cd ~/Soft/geant4.10.01.p03-install/bin/
    source geant4.sh
    cd -

#### ROOT installation

    sudo apt-get install git dpkg-dev cmake g++ gcc binutils libx11-dev libxpm-dev \
    libxft-dev libxext-dev

Also install optionally:

    sudo apt-get install gfortran libssl-dev libpcre3-dev \
    xlibmesa-glu-dev libglew1.5-dev libftgl-dev \
    libmysqlclient-dev libfftw3-dev libcfitsio-dev \
    graphviz-dev libavahi-compat-libdnssd-dev \
    libldap2-dev python-dev libxml2-dev libkrb5-dev \
    libgsl0-dev libqt4-dev

Download and compile ROOT:

    wget https://root.cern.ch/download/root_v6.10.08.source.tar.gz
    tar xf root_v6.10.08.source.tar.gz
    mkdir root_v6.10.08_build
    cd root_v6.10.08_build
    cmake ..
    make -j4

Add to your ~/.bashrc

    cd ~/Soft/root-6.10.08_build
    source bin/thisroot.sh
    cd -

#### OMGsim
OMGsim needs libconfig++

    sudo apt install libconfig++8-dev

Install SVN:

    sudo apt install subversion

Download the code and compile:

    git clone https://github.com/vkulikovskiy/OMGsim.git
    mkdir build
    cd build
    cmake ..
    make -j4

Then the binaries from the compilation will be in the folder "build".
The "programs" folder (might) contains the typical use of the tool, the test one contains some experimental programs.

##Quick start
A quick start howto is given [here](md_documentation_pages_quickstart.html).

##Documentation
###Build the documentation
The doxygen documentation is provided in the folder documentation.
It is easier to use it at your local machine.

Compile with:

    cd [your checkout path]/trunk
    doxygen documentation/config.dox
    cd documentation/latex
    make

Open (on Linux):

    firefox documentation/html/index.html
    gnome-open latex/refman.pdf

Open (on Mac):

    open documentation/html/index.html 
    open latex/refman.pdf

Note that the Doxygen needs to be installed. 

On Mac type:

    brew install doxygen  

###Complete the documentation
To complete the documentation you have to know some doxygen rules while you are comenting the .h code 
described [here](http://www.stack.nl/~dimitri/doxygen/manual/docblocks.html).

The easiest way could be to use the [doxymacs](https://github.com/emacsattic/doxymacs) script for emacs.
