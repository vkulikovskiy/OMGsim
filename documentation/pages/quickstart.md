Quick start
-----------

This short tutorial aims to introduce the using of OMGSim with some few examples.


##Make it run! now!

The base installation contains basics script and material/OM/PMT/detector descriptions and can be run immediatly.
The binary is OMGsim:

    ./OMGsim [-s] [-fa] ../script/vis-K40-DOM_local.mac [outputfile.root]

The -s is to open a GEANT4 session, unavoidable to draw the simulation (to be specified with your favorite way in the .mac file). The -fa is to get a "full analysis", just meaning that the output file will contain the full detected events (vertex, energy, impact point...) instead of the simplified events summary only (typicaly number of detected events)

##Custumize it a bit


###The detector data files


The vis-K40-DOM_local.mac contains the basics to do a K40 test with a DOM in the sea water. 
The detector data files are placed in common/data, should start with KM3Det, dans finish with .dat.
You can add any predefined volume (see further), DOM and PMT, the nomenclature should be like DOM1 DOM2 etc.

###The Passive Volumes

The passive volumes are determined by basics Geant4 solids or boolean operation. The used syntax is JSON.
Those files contain the name, the material and the volume kind.

###The DOM and OM
To build a DOM or OM, you have a bunch of parameter that allow you to define the size of the components and the materials.
[The parameters to define a DOM are listed here](structkm3net_1_1KM3OpticalModuleMgr_1_1OpticalModuleParameters.html) and should be implemented in JSON.
In the case of the DOM, the position and number of PMTs is not predifined, and is up to the PMTPositions table.

###The PMT

[The parameters to define a phototube](structkm3net_1_1KM3PhotoTubeMgr_1_1PhotoTubeParameters.html) usage is illustrated in the [3 inches](KM3PMT3inches_8dat_source.html).
![Caption text](../images/pmpar.svg "Image title")

###The Materials

Due to the "tabular like" structure of the material parameters an different configuration file structure is used.
The first element is the

    CREATE [name]
    [...]
    CREATE

Inside this the components flags, the density and components list are defined

    DENSITY [value in g/cm3]

    COMPONENTS
    [G4_component] relative density
    [...]
    COMPONENTS

Then the list of properties can be given if the form

    PROPERTY [name]
    [energy 1] [value 1]
    [energy 2] [value 2]
    ...

[The material of the photocathode](KM3MatPMT3inches_8dat_source.html) is a good example.
