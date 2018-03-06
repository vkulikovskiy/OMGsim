## Visualization
Visualization in Geant4 is still rather complicated. We suggest to use JAS3 for the visualization.
##JAS3
###Producing HepRep file on a machine running OMGsim.
First you will need a mac file with the HepRep visualization block (place it before /run/beamOn part):

    /vis/open HepRepXML
    /vis/scene/create filename.bheprep.zip
    /vis/drawVolume
    /vis/scene/endOfEventAction accumulate 500
    /tracking/storeTrajectory 1

Second Geant4 struggles with complicated design visualizations. Errors like this can be seen in the KM3NeT
DOM is visualized:

    ERROR: G4VSceneHandler::RequestPrimitives
    Polyhedron not available for 3inchesNoticeablePhotoCathode.
    This means it cannot be visualized on most systems.
    Contact the Visualization Coordinator.

It is suggested to reduce number of PMTs to one by using the "DOMAloneTest" detector and add the following
line somewhere after the detector name definition:

    /KM3/det/DontDraw bento absorber

The following complete example should work:

    /KM3/XP/Type K40
    /KM3/det/DetectorTypeName DOMAloneTest #K40cfg #DOMAlone

    /KM3/det/setTargetScattering 1
    /KM3/det/setTargetMaterial AntaresWater
    /KM3/det/DontDraw bento absorber
    /KM3/det/setWorldLength 10 m

    /KM3/GenType 1

    /run/verbose 0
    /run/initialize
    /tracking/storeTrajectory 0

    /gps/pos/shape Sphere
    /gps/pos/radius 10 m
    /gps/pos/centre 0 0 0
    /gps/pos/type Volume
    /gps/particle ion
    /gps/ion 19 40
    /gps/pos/confine Target

    /vis/open HepRepXML
    /vis/scene/create filename.bheprep.zip
    /vis/drawVolume
    /vis/scene/endOfEventAction accumulate 500
    /tracking/storeTrajectory 1

    /tracking/verbose 0
    /event/verbose 0
    /run/beamOn 1

After running OMGsim with a mandatory option "-s" the file scene-0.heprep.zip will be produced (the name
does not change even if filename.bheprep.zip  was requested).

    OMGsim -s your_mac.mac

If you want later to navigate through each event separately you can run Geant4 event by event placing in the 
mac file:

    /run/beamOn 1
    /run/beamOn 1
    /run/beamOn 1
    ...

Or introducing control loops with an additional mac file having "/run/beamOn 1".

###Installing JAS3 client on your local machine.
Java Runtime Engine should be installed to run JAS3 (Java Analysis Studio).
On Mac:

    brew cask install java

ON Ubuntu:

    sudo apt-get install default-jre default-jdk

Download and run latest JAS3 studio from http://jas.freehep.org/jas3/ -> Download -> jas-assembly-XXX-distribution.zip

    wget http://java.freehep.org/maven2/org/freehep/jas-assembly/3.2.7/jas-assembly-3.2.7-distribution.zip
    unzip jas-assembly-3.2.7-distribution.zip
    cd jas-assembly-3.2.7
    ./jas3

Install WIRED4 and "Jas 3 HepRep Plugin" going to View->Plugin Manager...->Available->Visualization 
(Select them, click Install selected plugins, restart JAS3).

File -> New -> Wired4 View

File -> Open File -> scene-0.heprep.zip

Press Play.

You can rotate, zoom unzoom. Select/unselect different parts of the detector unchecking boxes in 
Types->G4GeometryType.

Measure distances with a ruler tool.

Select different tracks by clicking on them with a picket tool (close to the ruler).
