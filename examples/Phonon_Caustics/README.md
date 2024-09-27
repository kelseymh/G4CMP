# G4CMP Phonon Simulation
Program to Create Caustic Phonon Plots.

In order to utilize qubits as particle detectors, understanding energy dissipation in qubit substrate (sapphire) phonon kinematics is essential. These mechanisms are strongly associated with correlated errors in qubit chips, as observed in cosmic muon and gamma ray absorption events reported recently. We attempt to better understand energy dissipation in sapphire with phonon caustics and kinematics simulation. For this, we expand the capabilities of G4CMP to include phonon kinematics in sapphire. In this G4CMP program, you can obtain phonon caustic images for sapphire, though it will also work for other GaAs, LiF, and SiO2.<br> 

Note: for the moment, kinematics for those substrates is only for phonon transport, and electron-hole pair transportation is still under construction.<br> 


For those interested, caustic measurement experimental methodology:

The substrate crystal has cooled to a low temperature (only a few degrees above absolute zero temperature). On one face of the substrate, a metal film is deposited and on the opposite face, a superconductor material (sensor) is deposited. The phonons are produced by passing a short burst of current through a metal or by hitting the metal with a source of focused laser light. After the excitation, several non-spherical shells of longitudinal and transverse phonons are produced and detected on the superconductor detector, increasing its temperature. It is possible to select the phonon mode by triggering the arrival time on the detector, phonons travel at different velocities on the substrate.

The ballistic phonons emitted from a point source concentrate along certain directions of the crystal, this is called phonon focusing.<br> 

If we have a small detector and we want to see how those lines are concentrated in the space we have two options. Move the detector to different positions in the space and later plot the intensity (number of lines detected on the detector) obtained at every point of the space. The other option is to keep the detector fixed and move the source point. In both cases, we are going to obtain the same plots of the intensities as a function of the position. In an experiment is not possible to move the detector, therefore the point source is moved. In the case of the simulation, we can  move both the phonon source or the detector, also create a big detector area that covers the maximum possible directions. <br> 

In our simulation, we follow the approach of creating a big detector, because it is significantly more computationally efficient.. 


 The phonon simulation is simple in terms of geometry. The basic elements of the simulation are 

* A sapphire substrate of 4mm x 4mm x 4mm.
* A big bolometer sensor of 4mm x 4mm x 0.0001 mm (Sensitive part where the phonons are absorbed).
* The new config.txt file with all the parameters for the sapphire substrate.
* The Caustic.mac file to generate 40e+6 phonos with isotropic initial random momentum.
* The Caustic_plot.root program to plot the caustic pattern for Transverse Slow Phonons, Transverse Fast Phonons, or both together.
# Substrate Dimension to Match the Experimental Phonon Caustics
The phonon caustics depends on the angular dimension and the size of the subtrate in exeprimental. <br> 
In order to match the experimental angles with the simualtion we need to initialize the initial phonon at difference distance <br> from the top part of the bolometers.
The angle is differents for each material. You need to use the following dimensions and the initial position of the phonons to reproduce the experimental phonon caustisc.<br> 
* Sapphire 4mm x 4mm x 4mm ```console /gps/pos/centre 0.0 0.0 0.08 cm ``` , &theta; = 58.0
* GaAs 4mm x 4mm x 4mm ```console /gps/pos/centre 0.0 0.0 0.08 cm ``` , &theta; = 59.04
* LiF 4mm x 4mm x 4mm ```console /gps/pos/centre 0.0 0.0 -0.03 cm ``` , &theta; = 40.0
* CaF2 4mm x 4mm x 8mm ```console /gps/pos/centre 0.0 0.0 -0.06 cm ``` , &theta; = 23.0
* CaWO4 9mm x 9mm x 6mm ```console /gps/pos/centre 0.0 0.0 0.0 cm ``` , &theta; = 56.3

# Miller Orientations
The phonon caustic depends on the orientation of the crystal with the sensor. You can obtain different phonon caustics images by changing the crystal orientation. <br>
In the case of Sapphire 4 Miller indices are needed. The program 4_Miller_to_3_Miller.py converts 4 Miller indices to 3 Miller.<br>
For more information check this link   [4 Miller indices ](https://apps.dtic.mil/sti/trecms/pdf/AD1115835.pdf)<br>
Note.-  The program 4_Miller_to_3_Miller.py is designed to use only for sapphire and match the [experimental results](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.29.2190) .

# How to Obtain the Phonon Caustics Patterns
We follow the standard procedure established on Geant4 to compile and execute a program. 
```console
User@fun:~$ cmake ..
User@fun:~$ make install
[100%] Built target g4cmpPhononCaustics
```
If you do not have errors
```console
./g4cmpaSapphire Caustic.mac
```
You will have a txt file "phonon_hits.txt" with three columns: <br> 
The first column is the name of the particle Transverse Fast or Transverse Slow <br> 
The second and third column is the final position (X,Y) of the phonons on the bolometer <br> 
We have an analysis root macro "Caustics_Plots.C" for plotting the output phonon_hits.txt. <br>
To plot only transverse fast Phonon  run the following command
```console
 root 'Caustics_Plots.C("Fast")'
```
for Transverse Slow
```console
root 'Caustics_Plots.C("Slow")'
```
Both Phonons
```console
root 'Caustics_Plots.C("Both")'
```
Your plot will look like 
for 0 0 0 1 Miller Orientation .



# Phonon Caustic Plots
![Alt text of the image](https://github.com/Israel-Tanjiro/Sapphire_G4CMP/blob/main/Sapphire_Phonon.png)
We do not include the analysis for longitudinal phonons, as they do not concentrate too much along certain directions of the crystal.<br>



# Using New Crystal Structures in Your Code
The folder crystal maps includes the config.txt Files for other substrate materials.
You can reproduce the phonon caustics pattern using the same program with the following lines changed in the Caustic_PhononDetectorConstruction.cc file. (You need to specify the materials.)
```ruby
fSubstrate = new G4Material("fSubtrate", 3.98*g/cm3, 2);
fSubstrate->AddElement(nistManager->FindOrBuildElement("Al"), 2);
fSubstrate->AddElement(nistManager->FindOrBuildElement("O"), 3);

```
The important part is the FindOrBuildElement and the Density of the Material.
```ruby
fSubstrate = new G4Material("fSubstrate", 5.32*g/cm3, 2);
fSubstrate->AddElement(nistManager->FindOrBuildElement("Ga"), 1);
fSubstrate->AddElement(nistManager->FindOrBuildElement("As"), 1);
```
The other line of the code that you must change is
```ruby
  G4LatticeLogical* SubstrateLogical = LM->LoadLattice(fSubstrate, "Al2O3");
```
to 

```ruby
  G4LatticeLogical* SubstrateLogical = LM->LoadLattice(fSubstrate, "GaAs");
```


The next step is to compile and execute the program using the standard cmake and make commands as outlined in the main G4CMP Readme.md file.


here’s a cool plot you can make for GaAs!
[experimental results](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.55.95)
![Alt text of the image](https://github.com/Israel-Tanjiro/Sapphire_G4CMP/blob/Substrate_G4CMP/Phonon_GaAS_110.png) 




