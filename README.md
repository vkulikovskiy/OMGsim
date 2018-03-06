# OMGsim
Optical Modules Geant4 simulation.

Authors: Christophe Hugon - CPPM (CNRS) and Vladimir Kulikovskiy - INFN Genova.

Documentation 
========================
The doxygen documentation is provided in the folder documentation.
It is easier to use it at your local machine.
Compile with:

    git clone https://github.com/vkulikovskiy/OMGsim
    cd OMGsim
    doxygen documentation/config.dox
    cd documentation/latex
    make

Open (on linux):

    firefox documentation/html/index.html
    gnome-open latex/refman.pdf

Open (on Mac):

    open documentation/html/index.html 
    open latex/refman.pdf

Note that the Doxygen needs to be installed. On MacOS type:

    brew install doxygen 
