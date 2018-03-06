##KM3NeT Digital Optical modules (DOMs)

##DOM definition

### PMT description
File: KM3PMT3inches.dat

PMT shape is defined from the Hamamatsu R12199-02 scheme.

* The front part of the PMT is a sphere with <RadiusSphere> internal radius of the PMT sphere (external is extended by the glass width). Chosen 50 + 2.5 mm.
* The thickest part of the PMT after the spherical part is approximated with an ellipsoid slice. <BulbRadiusMin>  = 35.95 mm is the radius of the ellispsoid slice adjusted to the sphericel part. The back side of the ellipsoid approximated with the same radius.
* <BulbRadiusMax>  40 mm is the thickest part of the ellipsoid/PMT.
* <BulbHeight>     18 mm is the thickness of the ellipsoid slice. It was obtained by measuring with a ruler Hamamatsu scheme "Dimension of R12199-02 Bulb-1-1.pdf
* <ConeHeight>     10 mm is the back conical shaped side of the PMT bulb.

### PMT collection efficiency calibration

The PMT collection efficiency is not simulated in current OMGsim version (and it is hardly possible it will be done in the future since it involves electron gain in dynode stack simulation). The collection efficiency was tuned according to the horizontal PMT scan performed in ECAP.

The simulation reproduces rather well the Quantum Efficiency measurements. The collection efficiency is stored in ANGULAR_EFFICIENCY field of common/data/KM3MatPMT3inches.dat as zenith angle on the PMT surface in rad vs efficiency. At first this was filled with efficency 1 everywhere and then the ratio between the obtained PMT efficiency in the simulation and Detection Efficiency in the simulation was  calculated using the routine below:

    ifstream ifs_scanx_theta("scanx_theta.txt");                                                          
    double angle_theta, scanx;                                                                            
    TGraph *gtheta2scanx = new TGraph();                                                                  
    while (ifs_scanx_theta >> scanx >> angle_theta)  gtheta2scanx->SetPoint(gtheta2scanx->GetN(),angle_theta,scanx);

    for (double angle = 0; angle <= 90; angle+=1) {                                                       
        double angle_theta = angle/180.*TMath::Pi()/3.*2.1;   //this is a bit weird convolution of the angles came from ANTARES and currently used to store the angle in ANGULAR_EFFICIENCY field 
        if (angle_theta<=gtheta2scanx->GetX()[gtheta2scanx->GetN()-1]) {                              
            double scanx = gtheta2scanx->Eval(angle_theta, 0, "S");                               
            double measDE;                                                                        
            if (scanx == 0) measDE = (0.191188+0.194693)/2.; //average by hand between two central points, value stored to 0 mm scan
            else measDE = (gDEmeas_pos->Eval(scanx, 0, "S")+gDEmeas_neg->Eval(scanx, 0, "S"))/2.;  //average between positive scan and negative scan
            //cout << scanx << " ";                                                               
            cout << angle_theta << " " << measDE/gAA2->Eval(scanx, 0, "S") << endl;  //the last sim graph is used here!
        }                                                                                             
    }

scanx_theta.txt stores the conversion between scan x (mm) and angle on the PMT. It was obtained with and output inside OMGsim code and commented out after.

    0 0
    10 0.14919
    20 0.305516
    25 0.389409
    30 0.479783
    35 0.581824
    37 0.636807
    38 0.683043

Due to the reflection on the back side of the PMT the obtained scan curve can stil mismatch the measurement. Angular Efficiency coefficients can be iteratively adjusted to match the data. 
### Light reflector ring
File: KM3OMkm3net.dat (or derivatives).

* The ring is obtained by rolling plate (flat ring sector with radial width of 18.05 cm). This width defines the area of the ring, so collected light and should be preserved.
* <ReflectorThickness> is the rolled reflector ring thickness (height of the cut cone shape) and is taken from the scheme (12 cm). Note, this is not the width of the unrolled ring!
* The opening angle is set to 48.33 degrees (according to private communication with Ronald Brujin the suggested value is 48 deg, however ReflectorThickness=12 cm can be obtained by this angle and 18.05 plate width). 
* The ring vertical width 13.5 cm from the scheme fits well 18.05 * cos(48.33) = 13.48.

Note that the reflector ring spacing from the glass is defined as a minimal distance (<ReflectorDistGlass> = 1.7 cm). For placement in Geant4 this distance is recalculated as a horisontal distance (assuming PMT is looking in the horizontal direction).

### DOM creation

The DOM creation is described in GetDOM function of common/src/KM3OpticalModuleMgr.cc. Note that PMTs are placed in polar coordinates (theta, phi) directrly in this function according to PMTPositions entry in common/data/KM3OMkm3net.dat. Their type and sizes are also read from the same file. In future release this may be done in different way (using \"Place them\" method in common/src/KM3DetectorConstruction.cc.

##K40 simulations

To run the K40 simulation run:

    cd build/scripts
    ../OMGsim -fa vis-K40-DOM_local.mac

And run the analysis:

    cd SimplePostAnalysis
    ./make.sh
    ./NoiseCoincidences ../build/prod/K40_DOMhadd.root coinc.root
    root -l calcK40rates.cc
