# OMGsim
Optical Modules Geant4 simulation

Documentation 
========================
The doxygen documentation is provided in the folder documentation.
It is easier to use it at your local machine.
Compile with:

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
