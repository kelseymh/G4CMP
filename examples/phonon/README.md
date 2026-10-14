# G4CMP Phonon propagation

> Daniel Brandt - SLAC (dbrandt@slac.stanford.edu)


This example demonstrates how phonon propagation in cryogenic crystals
can be simulated in G4CMP.

## Introduction

Phonon propagation is different from most other Geant4 propagation
simulations in a number of respects:

- Phonons are massless particles moving slower than the speed of light
- Phonon propagation and momentum vectors are not parallel
- Events isotropic in phonon-momentum space are not isotropic in real
 space.

This example will simulate the propagation of acoustic phonons through
a Germanium crystal, providing processes to simulate phonon scattering
of isotopic impurities, mode mixing between polarization states and
anharmonic downconversion (phonon splitting). As such it provides all
the physics required to realistically simulate phonon propagation in
cryogenically cold semiconductor crystals.

## Geometry

In this example the geometry is a cylindrical Germanium crystal
centered at (0,0,0) with Almuninium end caps. Phonons absorbed in the
Al end caps are counted by the sensitive detector.

## Primary event

The primary event is a single phonon of energy 7.5 meV at the center of
the Ge crystal. The polarization type (fast transvere, slow transverse or
longitudinal) is determined randomly according to the density of states
in Germanium. The direction of propagation is chosen randomly by the
Primary Generator Action class PhononPrimaryGeneratorAction.

## Execution and output

Source the G4CMP environment script, which sets the G4LATTICEDATA environment
variable so the crystal data files can be found.

```
source </path/to/your/G4CMP/install/share/G4CMP/g4cmp_env.sh>     # For SH/BASH users
source </path/to/your/G4CMP/install/share/G4CMP/g4cmp_env.csh>    # For CSH/TCSH users
```
Run the executable with one of the following commands.
```
./g4cmpPhonon                                                     # within the install directory
</path/to/your/Phonon_install/bin/g4cmpPhonon>                    # outside of the install directory
```

Upon execution, the usual Geant4 Qt interface opens. Run the macro `vis.mac` to
get a visualization of the Germanium detector. For the visualization to work,
OpenGL support must be installed.
The macro will automatically generate a single Primary Event (7.5 meV phonon)
at the center of the crystal.

The trajectory colour will indicate the polarization state of the phonon:

- Longitudinal (phononL):    blue
- Fast Transverse (phononTF): green
- Slow Transverse (phononTS): red

A small circle will be drawn wherever a phonon is absorbed into the 
Aluminium. All hits within the Aluminium are written into a
comma-separated-value (csv) file `phonon_hits.csv`, which can be opened with any
Calculator program (e.g. Excel, Libre office Calc) or with a plain text editor.

For these hits, the start energy, coordinates and time of the track is stored.
Also, the energy deposited into the Aluminium in the last step of the track is
listed together with the final coordinates and time of the post step point.

Every time a phonon is simulated, the information is appended to
`phonon_hits.csv`. If the file does not exist it will be created.

When running the `vis.mac`, screenshots of top and side of the Germanium crystal
are automatically saved in `g4cmpPhonon_0000.eps` and `g4cmpPhonon_0001.eps`.

The macro `single.mac` can be run either in the Geant4 Qt interface after
`vis.mac` to visualize 100 phonons (this might take a moment) or in batch mode
without visualization:

```
./g4cmpPhonon single.mac
```

Note that in `single.mac` the processes for phonon scattering (phononScattering)
and downconversion (phononDownconversion) are inactivated by default and can be
activated by commenting the corresponding macro commands by preceding the number
sign `#`. Be aware that activating phonon scattering increases the simulation
runtime significantly and a visualization may demand a lot of computing
resources.

In `single.mac` also the macro command `/tracking/verbose 1` is set, which lists
each step of each particle. This can produce a lot of output in the terminal.

## Visualization of an animation

If you would like to see an animation of the detector turning by 360 degrees,
open the Geant4 Qt interface, run `vis.mac` and afterwards `play.mac`. The
latter calls `loop.mac` which changes the viewing angle `theta` by steps of one
degree.
