# G4CMP -- Geant4 add-on framework for phonon and charge-carrier physics

    R. Agnese, D. Brandt, M. Kelsey, P. Redl


This package provide a collection of particle types, physics processes, and
supporting utilities to simulate a limited set of solid-state physics
processes in Geant4.  Developed for the low-temperature community, the
package support production and propagation of acoustic phonons and
electron-hole pairs through solid crystals such as germanium.

## Software Licenses

This product includes software developed by Members of the Geant4 
Collaboration ( http://cern.ch/geant4 ). A copy of the Geant4 license is
included at `G4CMP/Geant4\_LICENSE.html`.

This product ships an unaltered copy of Qhull source code. This is so that
we may compile and link the code in an unsupported way. The Qhull license
is included at `G4CMP/qhull-2012.1/COPYING.txt`.

The original parts of the product are licensed under the GNU General Public
License version 3 or (at your discretion) any later version. The full
license can be found at `G4CMP/LICENSE`.

## User Environment

Users must have a recent (10.2 or later) version of GEANT4 installed and
configured (via GEANT4's `bin/geant4.sh` or `bin/geant4.csh`. See GEANT4's
documentation for further instructions.).

Add the G4CMP environment variables using the `g4cmp\_env.csh` or `...sh`
scripts found in the G4CMP installation directory (see below for build and
installation procedures):

	source g4cmp\_env.csh		# For CSH/TCSH users
	. g4cmp\_env.sh			# For SH/BASH users

This must be done before building or running executables.

G4CMP is only configured for use on Linux and MacOSX platforms.  A minimum
configuration requires a recent enough version of GCC or Clang to support
the C++11 standard.

Several configuration parameters are available through environment variables
and macro commands, as listed below.  Most of these affect charge carrier
propagation and related processes. Because G4CMP is in a very active 
development state, do not expect this table to stay up-to-date. Rather,
developers should check the source code in 
`G4CMP/library/src/G4CMPConfigManager.cc` and 
`G4CMP/library/src/G4CMPConfigMessenger.cc` to see what is available.


| Environment variable    | Macro command                 | Value/action                            |
| ------------------------| ----------------------------- | ----------------------------------------|
| G4LATTICEDATA           | /g4cmp/LatticeData	          | Directory with lattice configs          |
| G4CMP\_DEBUG	           | /g4cmp/verbose <L> >0:        | Enable diagnostic messages              |
| G4CMP\_CLEARANCE [L]     | /g4cmp/clearance <L> mm       | Minimum distance of tracks from boundaries |
| G4CMP\_VOLTAGE [V]       | /g4cmp/voltage <V>	volt !=0:  | Apply uniform +Z voltage                |
| G4CMP\_EPOT\_FILE [F]     | /g4cmp/EPotFile <F> V=0:      | Read mesh field file "F"                |
| G4CMP\_EPOT\_SCALE [F]    | /g4cmp/scaleEPot <M> V=0:     | Scale the potentials in EPotFile by factor m|
| G4CMP\_MIN\_STEP [S]      | /g4cmp/minimumStep <S> S>0:   | Force minimum step S\*L0                |
| G4CMP\_MAKE\_PHONONS [R]  | /g4cmp/producePhonons <R>     | Fraction of phonons from energy deposit   |
| G4CMP\_MAKE\_CHARGES [R]  | /g4cmp/produceCharges <R>     | Fraction of charge pairs from energy deposit |
| G4CMP\_LUKE\_SAMPLE [R]   | /g4cmp/sampleLuke <R>         | Fraction of generated Luke phonons |
| G4CMP\_DOWN\_SAMPLE [R]   | /g4cmp/downconvertPhonons <R> | Fraction of natural phonon downconversions |
| G4CMP\_EMIN\_PHONONS [E]  | /g4cmp/minEPhonons <E> eV     | Minimum energy to track phonons         |
| G4CMP\_EMIN\_CHARGES [E]  | /g4cmp/minECharges <E> eV     | Minimum energy to track charges         |
| G4CMP\_USE\_KVSOLVER      | /g4mcp/useKVsolver [t|f]      | Use eigensolver for K-Vg mapping        |
| G4CMP\_MILLER\_H          | /g4cmp/orientation h k l      | Miller indices for lattice orientation  |
| G4CMP\_MILLER\_K          |                               |                                         |
| G4CMP\_MILLER\_L          |                               |                                         |
| G4CMP\_EH\_BOUNCES [N]    | /g4cmp/chargeBounces          | Maximum e/h reflections                 |
| G4CMP\_PHON\_BOUNCES [N]  | /g4cmp/phononBounces          | Maximum phonon reflections              |
| G4CMP\_HIT\_FILE [F]	    | /g4cmp/HitsFile <F>           | Write e/h hit locations to "F"          |

The default lattice orientation is to be aligned with the associated
G4VSolid coordinate system.  A different orientation can be specified by
setting the Miller indices (hkl) with `$G4CMP\_MILLER\_H`, `\_K`, and
`\_L`.

The environment variable `$G4CMP\_MAKE\_CHARGES` controls the rate (R) as a
fraction of total interactions, at which electron-hole pairs are produced
by energy partitioning.  Secondaries will be
produced with a track weight set to 1/R:

	unsetenv G4CMP_MAKE_CHARGES     # No new charge pairs generated
	setenv G4CMP_MAKE_CHARGES 1     # Generate e/h pair at every occurrence
	setenv G4CMP_MAKE_CHARGES 0.001 # Generate e/h pair 1:1000 occurrences

When secondary phonons are not produced, the equivalent energy is recorded as
non-ionizing energy loss (NIEL) on the track.  Generating seconary phonons
will significantly slow down the simulation.

The environment variable `$G4CMP\_MAKE\_PHONONS` controls the rate (R) as a
fraction of total interactions, at which "primary" phonons are produced (by
energy partitioning or recombination).  Secondaries will be produced with a
track weight set to 1/R:

	unsetenv G4CMP_MAKE_PHONONS     # No secondary phonons generated
	setenv G4CMP_MAKE_PHONONS 1     # Generate phonon at every occurrence
	setenv G4CMP_MAKE_PHONONS 0.001 # Generate phonon 1:1000 occurrences

When primary phonons are not produced, the equivalent energy is recorded as
non-ionizing energy loss (NIEL) on the track.  

Secondary phonons may be produced either by downconversion of higher energy
phonons, or by emission of Luke-Neganov phonons from charge carriers.
Generating seconary phonons can significantly slow down the simulation, so
each of these processes has analogous environment variables,
`$G4CMP\_LUKE\_SAMPLE` and `$G4CMP\_DOWN\_SAMPLE`, respectively, defined
with rates (R) as above.

For phonon propagation, a set of lookup tables to convert wavevector (phase
velocity) direction to group velocity are provided in the lattice
configuration file (see below).  The environment variable
`$G4CMP\_USE\_KVSOLVER` controls whether the eigenvalue solver should be
used directly for these calculations, instead of the lookup tables.  The
eigensolver imposes a factor of three penalty in CPU time, with the benefit
of maximum accuracy in phonon kinematics.

Three optional environment variables are used to configure the electric
field across the germanium crystal.  `$G4CMP\_VOLTAGE` specifies the voltage
across the crystal, used to generate a uniform electric field (no edge or
corner effects) from the bottom to the top face.  If the voltage is zero
(the default), then `$G4CMP\_EPOT\_FILE` specifies the name of the mesh
electric field field to be loaded for the g4cmpCharge test job.  There is no
default file.

For developers, there is a preprocessor flag (`make G4CMP\_DEBUG=1`) which may
be set before building the libraries.  This variable will turn on some
additional diagnostic output files which may be of interest.

## Building the Package

G4CMP supports building itself with either GNU Make or CMake, and separately
supports being linked into user applicated with either GNU Make (via
environment variable settings) or CMake.

### Building with Make

Configure your Geant4 build environment using
`<g4dir>/share/Geant4-${VERSION}/geant4make/geant4make.csh` or `...sh`, then
configure or G4CMP environment as described above with `g4cmp_env.csh` or
`...sh`.

After configuring your environment, build the G4CMP library with the command

	make library

The libraries (libg4cmp.so and libqhull.so) will be written to your
`$G4WORKDIR/lib/$G4SYSTEM/` directory, just like any other Geant4 example or
user code, and should be found automatically when linking an application.

With the library built, any of the three demonstration programs (phonon,
charge) may be built as a normal GEANT4 user application.
From the top-level directory, use the command

	make examples

to build them all, or

	make <name>

to build just one (where <name> is the directory name of interest).  The
executables will be named `g4cmpPhonon` and `g4cmpCharge`, respectively, and
will be written to `$G4WORKDIR/bin/$G4SYSTEM/`.

### Building with CMake

Create a build directory outside of the source tree,

    mkdir /path/to/G4CMP/../G4CMP-build
    cd /path/to/G4CMP-build

We must tell CMake where GEANT4 is installed. If you want only the library
to be built, use the following command

    cmake -DGeant4_DIR=/path/to/Geant4/lib64/Geant4-${VERSION} ../G4CMP

If you want to install to a local path, rather than system-wide, use the
`-DCMAKE\_INSTALL\_PREFIX=/path/to/install` option.

If you want to build an example application,

    cmake -DGeant4_DIR=/path/to/Geant4/lib64/Geant4-${VERSION} -DBUILD_CHARGE_EXAMPLE=ON ../G4CMP

If you want to build all examples,

    cmake -DGeant4_DIR=/path/to/Geant4/lib64/Geant4-${VERSION} -DBUILD_ALL_EXAMPLES=ON ../G4CMP

Then simply run the `make` command in the build directory

    make

While it's not strictly necessary, we strongly recommend installing G4CMP to 
the install prefix rather than running the binaries from the build directory

    make install

Once theinstall step is completed, the /path/to/install/share/G4CMP/
directory will contain copies of the `g4cmp\_env.csh` and `...sh` scripts
discussed above.  These copies should be sourced in order to correctly
locate the installed libraries and header files.

### Linking user applications against G4CMP

G4CMP is an application library, which can be linked into a user's Geant4
application in order to provide phonon and charge carrier transport in
crystals.  Users must reference G4CMP in their application build in order to
utilize these features.  

If you have a simple Makefile build system (GMake), the following two lines,
or an appropriate variation on them, should be sufficient:

    CXXFLAGS += -I$(G4CMPINCLUDE)
    LDFLAGS += -L$(G4CMPLIB) -lG4cmp -lqhullcpp -lqhullstatic_p

These actions must occur _before_ the Geant4 libraries and include directory
are referenced (G4CMP includes modified versions of some toolkit code).

If you are using CMake to build your application, it should be sufficient to
add the following two actions, before referencing Geant4:

    find_package(G4CMP REQUIRED)
    include(${G4CMP_USE_FILE})


## Defining the Crystal Dynamics

In a user's application, each active G4CMP material (e.g., germanium
crystals or diamonds) must have a collection of dynamical parameters
defined.  These parameters are used by the phonon and charge-carrier
processes to know how to create, propagate, and scatter the particles
through the crystal.

Each material's parameters are stored in a subdirectory under CrystalMaps
(or wherever the envrionment variable `$G4LATTICEDATA` points).  G4CMP is
distributed with germanium and silicon configurations, in `CrystalMaps/Ge/`
and `CrystalMaps/Si/`, respectively.  We recommend naming additional
directories by element or material, matching the Geant4 conventions, but
this is not required or enforced.

The parameter definition file is config.txt.  Each line starts with a
keyword, followed by one or more values.  Any line which starts with "#" is
ignored, as is any text after a "#" on a parameter line.  Multiple
keywords/value sets may be included on a single line of the file if desired
for readability.

Dimensional parameters MUST be specified with the value in each entry.  For 
keywords taking multiple values, a single unit may be specified after the
group of values, e.g., 

	triclinic 1. 2. 3. Ang 30. 50. 45. deg

where "Ang" and "deg" are the appropriate length and angular dimensions.

The lattice symmetry is specified by one of the seven crystal systems (or
"amorphous") followed by the appropriate combination of lattice constant(s)
and angle(s) needed to specify it uniquely.  The reduced elasticity matrix,
Cij, must be specified term by term; which components are needed depends on
the crystal system.

| Keyword | Arguments | Value type(s)             | Units              |
|---------|-----------|---------------------------|--------------------|
| **Lattice parameters** |
| amorphous | -none-  | Polycrystalline solid     |                    |
| cubic   | a         | Lattice constant          | length             |
| tetragonal | a c    | Lattice constants         | length             |
| hexagonal  | a c    | Lattice constants         | length             |
| orthorhombic | a b c | Lattice constants        | length             |
| rhombohedral | a alpha | Lattice const., angle  | length, deg/rad    |
| monoclinic | a b c alpha | Lattice const., angle | length, deg/rad   |
| triclinic | a b c alpha beta gamma | Lattice const., angle | length, deg/rad |
| stiffness | i j val | Indices 1-6, elasticity   | pressure (Pa, GPa) |
| Cij       | i j val | Indices 1-6, elasticity   | Pa, GPa            |
| **Phonon parameters** |
| beta    | val       | scattering parameters     | Pa, GPa            |
| gamma   | val       | (see S. Tamura, PRB 1985) | Pa, GPa            |
| lambda  | val       |                           | Pa, GPa            |
| mu      | val       |                           | Pa, GPa            |
| dyn    | beta gamma lambda mu | All four params | Pa, GPa            |
| scat    | B         | isotope scattering rate   | second^3 (s3)      |
| decay   | A         | anharmonic decay rate     | second^4 (s4)      |
| decayTT | frac      | Fraction of L->TT decays  |                    |
| LDOS    | frac      | longitudnal density of states | sum to unity   |
| STDOS   | frac      | slow-transverse density of states |            |
| FTDOS   | frac      | fast-transverse density of states |            |
| Debye   | val       | Debye energy for phonon primaries | E, T, Hz   |
| **Charge carrier parameters** |
| vsound  | Vlong     | sound speed (longitudinal) | m/s               |
| vtrans  | Vtrans    | sound speed (transverse)   | m/s               |
| bandgap | val       | Bandgap energy             | energy (eV)       |
| pairEnergy | val    | Energy taken by e-h pair   | energy (eV)       |
| fanoFactor | val    | Spread of e-h pair energy  |                   |
| l0_e    | len       | electron scattering length | length            |
| l0_h	  | len       | hole scattering length     | length            |
| hmass   | m_h       | effective mass of hole   | electron mass ratio |
| emass   | m_xx m_yy m_zz | electron mass tensor | (same)             |
| valley  | theta phi psi unit | Euler angles     | angle (deg/rad)    |
| ivField | E0        | Edelweiss intervalley scattering  | V/m        |
| ivRate  | freq      | Intervalley scattering rate       | Hz         |
| ivPower | val       | Power law [sqrt(field^ivPower)]   | none       |
