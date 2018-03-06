##Calibration

##QE calibration

The most direct way to measure QE is with Direct Coupling to the PMT photocathode (DC). 
In this case every escaped electron participates in the signal.
The similar way is to measure the current from the first dynode. KM3NeT ECAP group reports negligible difference (99% vs 100%) for the two method comparison for the 3-inch PMTs. The illumination of the PMT should be done to the central part with rather small beam size (~1 mm) to minimise reflection and other effects for the inclinated beams falling on the spherical PMT surface.

Alternative is the measure of the QE from the amplified anode signal. In this case the collection efficiency (CE) is also included in the measurement. To calculate the QE some collection efficiency should be assumed (90% is suggested for big PMTs).

In case of ANTARES PMTs the anode signal was used (QE+CE). In case of KM3NeT PMTs the DC is used (QE).

Starting from v2.0 the simulation uses complex refractive index (N=n+ik). Known n, k for every wavelength and the photocathode thickness allows absorption, transmission and reflection calculation for every wavelength and the inclination angle. The "k" index can be deduced from QE measurement if the photocathode thickness is known. 

The real thickness is not provided officially by Hamamatsu. Possibilities:
- 30 nm 3-inch, 26 nm >= 8-inch (unofficial)
- 20 nm best fit of Motta et al.
- 26 nm typical of M.D.Lay meas. For 8-inch?

It was also verified that for every photocathode thickness (20,26,30) the measured average QE for 3-inch KM3NeT PMTs can be tuned. For the following it was decided to use 30 nm thickness for 3 inch PMTs.

To run the calibration, compile OMGsim with ENABLE_TESTING:

    cd build
    cmake .. -DENABLE_TESTING=1
    make

Run QuantumEfficiencyTester. Without options, it reads curring K-index in the data, simulates QE with decired precision (5% is the default) and compares it with the measured one. The precision comes from the electron production probability, which is a random Brownian process (this probability depends on the photocathod thickness and light penetration probability which is the function of k-index and wavelength, so it needs to be simulated with many electrons for every k-index). Precision can be changed with "-precision option". To calibrate the K-index, "-cal" option should be used. After several itarations (untill ratio of simulated QE to the measured stays inside the requested precision for every wavelength) the K-index table is printed. The following command was used for 3-inch KM3NeT calibration with 1% precision:

    ./QuantumEfficiencyTester -precision 0 -cal

### KM3NeT 3-inch PMT calibration

The result of the K-index calibration for Hamamatsu 3-inche PMTs (R12199-02) is shown below. Note, that the PMT glass was adjusted in a way that in UV part (280-350 nm) the transmission is almost 100%. 

    current KINDEX
                                                                                                                      2.1683           |0
           250 |0
           260 |0
           270 |0
           280 |0
           290 |---0.061594
           300 |---------0.18549
           310 |--------------------0.42085
           320 |----------------------------------0.73516
           330 |----------------------------------------------------1.1286
           340 |------------------------------------------------------------------1.4152
           350 |------------------------------------------------------------------------1.5541
           360 |----------------------------------------------------------------------------1.6457
           370 |---------------------------------------------------------------------------------------1.8901
           380 |----------------------------------------------------------------------------------------------------2.1683
           390 |----------------------------------------------------------------------------------------------------2.1654
           400 |-----------------------------------------------------------------------------------------------2.0531
           410 |---------------------------------------------------------------------------------------1.881
           420 |---------------------------------------------------------------------------------1.7474
           430 |-----------------------------------------------------------------------------1.6717
           440 |------------------------------------------------------------------------1.5612
           450 |----------------------------------------------------------------1.3731
           460 |-------------------------------------------------------1.1945
           470 |---------------------------------------------0.97649
           480 |---------------------------------------0.83466
           490 |----------------------------------0.73765
           500 |-------------------------------0.67252
           510 |----------------------------0.60474
           520 |------------------------0.51648
           530 |------------------0.38568
           540 |-------------0.27831
           550 |----------0.21523
           560 |---------0.17847
           570 |--------0.15562
           580 |-------0.13702
           590 |------0.12103
           600 |-----0.10712
           610 |-----0.095318
           620 |----0.083897
           630 |----0.073658
           640 |---0.063215
           650 |---0.054701
           660 |---0.049425
           670 |--0.039306
           680 |--0.029667
           690 |-0.0208
           700 |-0.01852
           710 |-0.011035
           720 |-0.0038807
           730 |0
           740 |0
    QE calculated from current KINDEX
                                                                                                                     0.27013           |0
           250 |0
           260 |0
           270 |0
           280 |-4.3213e-08
           290 |------------0.029886
           300 |-----------------------------0.077076
           310 |---------------------------------------------------0.13633
           320 |---------------------------------------------------------------------0.18458
           330 |----------------------------------------------------------------------------------0.22038
           340 |-----------------------------------------------------------------------------------------0.24022
           350 |--------------------------------------------------------------------------------------------0.24828
           360 |----------------------------------------------------------------------------------------------0.25187
           370 |-------------------------------------------------------------------------------------------------0.26005
           380 |----------------------------------------------------------------------------------------------------0.26956
           390 |----------------------------------------------------------------------------------------------------0.27013
           400 |---------------------------------------------------------------------------------------------------0.26676
           410 |-------------------------------------------------------------------------------------------------0.26147
           420 |-----------------------------------------------------------------------------------------------0.25683
           430 |----------------------------------------------------------------------------------------------0.25201
           440 |--------------------------------------------------------------------------------------------0.24714
           450 |----------------------------------------------------------------------------------------0.23818
           460 |-------------------------------------------------------------------------------------0.22799
           470 |-------------------------------------------------------------------------------0.21131
           480 |-------------------------------------------------------------------------0.19744
           490 |---------------------------------------------------------------------0.18571
           500 |------------------------------------------------------------------0.17793
           510 |---------------------------------------------------------------0.16801
           520 |---------------------------------------------------------0.1537
           530 |------------------------------------------------0.12781
           540 |--------------------------------------0.10245
           550 |--------------------------------0.084453
           560 |----------------------------0.0734
           570 |-------------------------0.065517
           580 |----------------------0.058418
           590 |--------------------0.052383
           600 |------------------0.047088
           610 |----------------0.042121
           620 |--------------0.037404
           630 |-------------0.033051
           640 |-----------0.028541
           650 |----------0.024767
           660 |---------0.022424
           670 |-------0.017857
           680 |------0.013645
           690 |----0.0095716
           700 |----0.0084958
           710 |--0.0050802
           720 |-0.0017878
           730 |0
           740 |0
    QE in the data
                                                                                                                     0.26912           |0
           250 |0
           260 |0
           270 |0
           280 |0
           290 |------------0.029973
           300 |-----------------------------0.077265
           310 |---------------------------------------------------0.13637
           320 |--------------------------------------------------------------------0.18309
           330 |----------------------------------------------------------------------------------0.22088
           340 |-----------------------------------------------------------------------------------------0.23947
           350 |--------------------------------------------------------------------------------------------0.24691
           360 |----------------------------------------------------------------------------------------------0.25115
           370 |-------------------------------------------------------------------------------------------------0.26125
           380 |----------------------------------------------------------------------------------------------------0.26842
           390 |----------------------------------------------------------------------------------------------------0.26912
           400 |---------------------------------------------------------------------------------------------------0.26643
           410 |--------------------------------------------------------------------------------------------------0.26202
           420 |------------------------------------------------------------------------------------------------0.25784
           430 |----------------------------------------------------------------------------------------------0.25277
           440 |--------------------------------------------------------------------------------------------0.24629
           450 |-----------------------------------------------------------------------------------------0.23798
           460 |------------------------------------------------------------------------------------0.22622
           470 |------------------------------------------------------------------------------0.21025
           480 |-------------------------------------------------------------------------0.19657
           490 |---------------------------------------------------------------------0.18507
           500 |------------------------------------------------------------------0.17681
           510 |---------------------------------------------------------------0.16765
           520 |---------------------------------------------------------0.15302
           530 |------------------------------------------------0.12798
           540 |--------------------------------------0.10241
           550 |--------------------------------0.084572
           560 |----------------------------0.07289
           570 |-------------------------0.065143
           580 |----------------------0.05849
           590 |--------------------0.052391
           600 |------------------0.046967
           610 |----------------0.042119
           620 |--------------0.037433
           630 |-------------0.032911
           640 |-----------0.02853
           650 |----------0.024741
           660 |---------0.022358
           670 |-------0.017918
           680 |------0.013605
           690 |----0.0095899
           700 |----0.0085087
           710 |--0.0050991
           720 |-0.0017842
           730 |0
           740 |0
    wavelength kindex 
    250 0
    260 0
    270 0
    280 0
    290 0.061594
    300 0.18549
    310 0.42085
    320 0.73516
    330 1.1286
    340 1.4152
    350 1.5541
    360 1.6457
    370 1.8901
    380 2.1683
    390 2.1654
    400 2.0531
    410 1.881
    420 1.7474
    430 1.6717
    440 1.5612
    450 1.3731
    460 1.1945
    470 0.97649
    480 0.83466
    490 0.73765
    500 0.67252
    510 0.60474
    520 0.51648
    530 0.38568
    540 0.27831
    550 0.21523
    560 0.17847
    570 0.15562
    580 0.13702
    590 0.12103
    600 0.10712
    610 0.095318
    620 0.083897
    630 0.073658
    640 0.063215
    650 0.054701
    660 0.049425
    670 0.039306
    680 0.029667
    690 0.0208
    700 0.01852
    710 0.011035
    720 0.0038807
    730 0
    740 0
    Ratio of calculated QE to QE in the data
                                                                                                                      1.0081           |0
           250 |0
           260 |0
           270 |0
           280 |0
           290 |---------------------------------------------------------------------------------------------------0.9971
           300 |---------------------------------------------------------------------------------------------------0.99756
           310 |---------------------------------------------------------------------------------------------------0.99968
           320 |----------------------------------------------------------------------------------------------------1.0081
           330 |---------------------------------------------------------------------------------------------------0.99777
           340 |----------------------------------------------------------------------------------------------------1.0031
           350 |----------------------------------------------------------------------------------------------------1.0055
           360 |----------------------------------------------------------------------------------------------------1.0029
           370 |---------------------------------------------------------------------------------------------------0.99538
           380 |----------------------------------------------------------------------------------------------------1.0043
           390 |----------------------------------------------------------------------------------------------------1.0038
           400 |----------------------------------------------------------------------------------------------------1.0012
           410 |---------------------------------------------------------------------------------------------------0.99789
           420 |---------------------------------------------------------------------------------------------------0.99608
           430 |---------------------------------------------------------------------------------------------------0.99698
           440 |----------------------------------------------------------------------------------------------------1.0035
           450 |----------------------------------------------------------------------------------------------------1.0009
           460 |----------------------------------------------------------------------------------------------------1.0078
           470 |----------------------------------------------------------------------------------------------------1.005
           480 |----------------------------------------------------------------------------------------------------1.0044
           490 |----------------------------------------------------------------------------------------------------1.0035
           500 |----------------------------------------------------------------------------------------------------1.0063
           510 |----------------------------------------------------------------------------------------------------1.0022
           520 |----------------------------------------------------------------------------------------------------1.0044
           530 |---------------------------------------------------------------------------------------------------0.99873
           540 |----------------------------------------------------------------------------------------------------1.0004
           550 |---------------------------------------------------------------------------------------------------0.99859
           560 |----------------------------------------------------------------------------------------------------1.007
           570 |----------------------------------------------------------------------------------------------------1.0057
           580 |---------------------------------------------------------------------------------------------------0.99877
           590 |---------------------------------------------------------------------------------------------------0.99985
           600 |----------------------------------------------------------------------------------------------------1.0026
           610 |---------------------------------------------------------------------------------------------------1
           620 |---------------------------------------------------------------------------------------------------0.99922
           630 |----------------------------------------------------------------------------------------------------1.0042
           640 |----------------------------------------------------------------------------------------------------1.0004
           650 |----------------------------------------------------------------------------------------------------1.0011
           660 |----------------------------------------------------------------------------------------------------1.0029
           670 |---------------------------------------------------------------------------------------------------0.99663
           680 |----------------------------------------------------------------------------------------------------1.003
           690 |---------------------------------------------------------------------------------------------------0.9981
           700 |---------------------------------------------------------------------------------------------------0.99849
           710 |---------------------------------------------------------------------------------------------------0.99629
           720 |----------------------------------------------------------------------------------------------------1.002
           730 |0
           740 |0
    precision 0.01 reached

The provided table "wavelength kindex" was inserted in the common/data/KM3MatPMT3inches.dat file, section "PROPERTY  KINDEX".