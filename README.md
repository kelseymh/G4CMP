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
included at `G4CMP/Geant4_LICENSE.html`.

This product ships an unaltered copy of Qhull source code. This is so that
we may compile and link the code in an unsupported way. The Qhull license
is included at `G4CMP/qhull-2012.1/COPYING.txt`.

The original parts of the product are licensed under the GNU General Public
License version 3 or (at your discretion) any later version. The full
license can be found at `G4CMP/LICENSE`.

## User Environment

Users must have a recent (10.4 through 10.7) version of GEANT4 installed and
configured (via GEANT4's `bin/geant4.sh` or `bin/geant4.csh`. See GEANT4's
documentation for further instructions.).

**NOTE** The relase of Geant4 Version 11 introduced substantial and breaking
  changes to many Geant4 interface classes.  We are maintaining G4CMP under
  ==Geant4 Version 10== (through 10.7) to ensure compatibility with our
  major experimental users.

Add the G4CMP environment variables using the `g4cmp_env.csh` or `...sh`
scripts found in the G4CMP installation directory (see below for build and
installation procedures):

	source g4cmp_env.csh		# For CSH/TCSH users
	. g4cmp_env.sh			# For SH/BASH users

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
| G4CMP\_DEBUG	          | /g4cmp/verbose [L] >0:        | Enable diagnostic messages              |
| G4CMP\_CLEARANCE [L]    | /g4cmp/clearance [L] mm       | Minimum distance of tracks from boundaries |
| G4CMP\_VOLTAGE [V]      | /g4cmp/voltage [V]	volt !=0:  | Apply uniform +Z voltage                |
| G4CMP\_EPOT\_FILE [F]   | /g4cmp/EPotFile [F] V=0:      | Read mesh field file "F"                |
| G4CMP\_EPOT\_SCALE [F]  | /g4cmp/scaleEPot [M] V=0:     | Scale the potentials in EPotFile by factor m|
| G4CMP\_MIN\_STEP [S]    | /g4cmp/minimumStep [S] S>0:   | Force minimum step S\*L0                |
| G4CMP\_EH\_BOUNCES [N]    | /g4cmp/chargeBounces [N]      | Maximum e/h reflections                 |
| G4CMP\_PHON\_BOUNCES [N]  | /g4cmp/phononBounces [N]      | Maximum phonon reflections              |
| G4CMP\_MAKE\_PHONONS [R] | /g4cmp/producePhonons [R]     | Fraction of phonons from energy deposit   |
| G4CMP\_MAKE\_CHARGES [R] | /g4cmp/produceCharges [R]     | Fraction of charge pairs from energy deposit |
| G4CMP\_LUKE\_SAMPLE [R] | /g4cmp/sampleLuke [R]         | Fraction of generated Luke phonons |
| G4CMP\_MAX\_LUKE [N] | /g4cmp/maxLukePhonons [N] | Soft maximum Luke phonons per event |
| G4CMP\_SAMPLE\_ENERGY [E] | /g4cmp/samplingEnergy [E] eV  | Energy above which to downsample |
| G4CMP\_COMBINE\_STEPLEN [L] | /g4cmp/combiningStepLength [L] mm | Combine
hits below step length |
| G4CMP\_EMIN\_PHONONS [E] | /g4cmp/minEPhonons [E] eV     | Minimum energy to track phonons         |
| G4CMP\_EMIN\_CHARGES [E] | /g4cmp/minECharges [E] eV     | Minimum energy to track charges         |
| G4CMP\_USE\_KVSOLVER    | /g4mcp/useKVsolver [t\|f]     | Use eigensolver for K-Vg mapping        |
| G4CMP\_FANO\_ENABLED    | /g4cmp/enableFanoStatistics [t\|f] | Apply Fano statistics to input ionization |
| G4CMP\_IV\_RATE\_MODEL  | /g4cmp/IVRateModel [IVRate\|Linear\|Quadratic] | Select intervalley rate parametrization |
| G4CMP\_ETRAPPING\_MFP   | /g4cmp/eTrappingMFP [L] mm        | Mean free path for electron trapping |
| G4CMP\_HTRAPPING\_MFP   | /g4cmp/hTrappingMFP [L] mm        | Mean free path for charge hole trapping |
| G4CMP\_EDTRAPION\_MFP | /g4cmp/eDTrapIonizationMFP [L] mm | MFP for e-trap ionization by e- |
| G4CMP\_EATRAPION\_MFP | /g4cmp/eATrapIonizationMFP [L] mm | MFP for h-trap ionization by e- |
| G4CMP\_HDTRAPION\_MFP | /g4cmp/hDTrapIonizationMFP [L] mm | MFP for e-trap ionization by h+ |
| G4CMP\_HATRAPION\_MFP | /g4cmp/hATrapIonizationMFP [L] mm | MFP for h-trap ionization by h+ |
| G4CMP\_TEMPERATURE   | /g4cmp/temperature [T] K | Device/substrate/etc. temperature |
| G4CMP\_NIEL\_FUNCTION | /g4cmp/NIELPartition [LewinSmith\|Lindhard] | Select NIEL partitioning function |
| G4CMP\_CHARGE\_CLOUD     | /g4cmp/createChargeCloud [t\|f] | Create charges in sphere around location |
| G4CMP\_MILLER\_H          | /g4cmp/orientation [h] [k] [l] | Miller indices for lattice orientation  |
| G4CMP\_MILLER\_K          |                               |                                         |
| G4CMP\_MILLER\_L          |                               |                                         |
| G4CMP\_HIT\_FILE [F]	    | /g4cmp/HitsFile [F]           | Write e/h hit locations to "F"          |

The default lattice orientation is to be aligned with the associated
G4VSolid coordinate system.  A different orientation can be specified by
setting the Miller indices (hkl) with `$G4CMP_MILLER_H`, `_K`, and
`_L`.

The environment variable `$G4CMP_MAKE_CHARGES` controls the rate (R) as a
fraction of total interactions, at which electron-hole pairs are produced
by energy partitioning.  Secondaries will be
produced with a track weight set to 1/R:

	unsetenv G4CMP_MAKE_CHARGES     # No new charge pairs generated
	setenv G4CMP_MAKE_CHARGES 1     # Generate e/h pair at every occurrence
	setenv G4CMP_MAKE_CHARGES 0.001 # Generate e/h pair 1:1000 occurrences

When secondary phonons are not produced, the equivalent energy is recorded as
non-ionizing energy loss (NIEL) on the track.  Generating seconary phonons
will significantly slow down the simulation.

The environment variable `$G4CMP_MAKE_PHONONS` controls the rate (R) as a
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
the `LukeScattering` process has an analogous environment variable,
`$G4CMP_LUKE_SAMPLE`, defined with rate (R) as above.

For simulations which generate primary phonons and charge carriers from
Geant4 energy deposition (using `G4CMPEnergyPartition`), the above
environment variables may be replaced with a sampling "energy scale,"
`$G4CMP_SAMPLE_ENERGY`.  This parameter is applied to each energy deposit,
and to ionization or NIEL energy separately.  If the energy deposit is below
the scale, then no biasing will be done (the scale factors will all be set
to 1.).  Above the energy scale setting, the scale factors will be set
according to E_scale_/E_deposit_.  In this mode, an additional sampling
parameter, `$G4CMP_MAX_LUKE`, may be set.  This sets an approximate maximum
number of Luke-Neganov phonons to be produced per event; the default is
about 10,000.

The parameter `$G4CMP_COMBINE_STEPLEN` (`/g4cmp/combiningStepLength`)
specifies a minimum step length for individual `G4CMPEnergyPartition` hits.
Shorter contiguous steps by a track will be consolidated into one hit, which
will then be processed with sampling as described above.  When combined with
Geant4's built in secondary production cuts, this should improve runtime
performance substantially.

For phonon propagation, a set of lookup tables to convert wavevector (phase
velocity) direction to group velocity are provided in the lattice
configuration file (see below).  The environment variable
`$G4CMP_USE_KVSOLVER` controls whether the eigenvalue solver should be
used directly for these calculations, instead of the lookup tables.  The
eigensolver imposes a factor of three penalty in CPU time, with the benefit
of maximum accuracy in phonon kinematics.

Three optional environment variables are used to configure the electric
field across the germanium crystal.  `$G4CMP_VOLTAGE` specifies the voltage
across the crystal, used to generate a uniform electric field (no edge or
corner effects) from the bottom to the top face.  If the voltage is zero
(the default), then `$G4CMP_EPOT_FILE` specifies the name of the mesh
electric field field to be loaded for the g4cmpCharge test job.  There is no
default file.

For developers, there is a preprocessor flag (`make G4CMP_DEBUG=1`) which may
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

If you want debugging symbols included with the G4CMP library, you
need to build with the G4DEBUG environment or Make variable set:

	export G4DEBUG=1
or
	setenv G4DEBUG 1
or
	make library G4DEBUG=1

If you want to enable additional diagnostics in some processes, including
writing out statistics files, build with the G4CMP_DEBUG environment or Make
variable set.  Note that this is not compatible with running multiple worker
threads.

	export G4CMP_DEBUG=1
or
	setenv G4CMP_DEBUG 1
or
	make library G4CMP_DEBUG=1

If you want to enable "sanitizing" options with the library, to look for
memory leaks, thread collisions etc., you may set the options
G4CMP_USE_SANITIZER and G4CMP_SANITIZER_TYPE (default is "thread"):

	export G4CMP_USE_SANITIZER=1
or
	setenv G4CMP_USE_SANITIZER 1
or
	make library G4CMP_USE_SANITIZER=1

*NOTE*:  If your source directory was not cloned from GitHub (specifically,
if it does not contain `.git/`) you may need to specify a version string for
identify the G4CMP version at runtime.  Use `G4CMP_VERSION=X.Y.Z` on the
Make command line for this purpose.  If `.git/` is available, the option
will be ignored.

### Building with CMake

Create a build directory outside of the source tree, such as

    mkdir /path/to/G4CMP/../G4CMP-build
    cd /path/to/G4CMP-build

We must tell CMake where GEANT4 is installed. If you want only the library
to be built, use the following command

    cmake -DGeant4_DIR=/path/to/Geant4/lib64/Geant4-${VERSION} ../G4CMP

By default, CMake will install a software package under /usr/local.  If you
want to install to a local path, rather than system-wide, use the
`-DCMAKE_INSTALL_PREFIX=/path/to/install` option.

If you want debugging symbols included with the G4CMP library, you
need to include the `-DCMAKE_BUILD_TYPE=Debug` option.

If you want to enable additional diagnostics in some processes, including
writing out statistics files, include the `-DG4CMP_DEBUG=1` option.  Note
that this is not compatible with running multiple worker threads.

If you want to enable "sanitizing" options with the library, to look for
memory leaks, thread collisions etc., you may set the options
`-DG4CMP_USE_SANITIZER=ON` and (optionally) `-DG4CMP_SANITIZER_TYPE=value`
(default is "thread", other values may be "memory", "address", or "leak").
If you do this, we recommend using the "Debug" build type.

**NOTE**:  If your source directory was not cloned from GitHub (specifically,
if it does not contain `.git/`) you may need to specify a version string for
identify the G4CMP version at runtime.  Use the `-DG4CMP_VERSION=X.Y.Z`
option for this purpose.  If `.git/` is available, the option will be
ignored.

If you want to copy the examples directories (see below) to the installation
area, use the option `-DINSTALL_EXAMPLES=ON` (for all examples). Each example
has been set up as a standalone "project" for CMake and can be configured via:

    cmake -DGeant4_DIR=/path/to/Geant4/lib64/Geant4-${VERSION} -DINSTALL_EXAMPLES=ON ../G4CMP

Once you've configured the build with `cmake` and option flags, run the
`make` command in the build directory

    make

and transfer the successfully build libraries to your installation area

    make install

Once the install step is completed, the /path/to/install/share/G4CMP/
directory will contain copies of the `g4cmp_env.csh` and `...sh` scripts
discussed above.  These copies should be sourced in order to correctly
locate the installed libraries and header files.

### Linking user applications against G4CMP

G4CMP is an application library, which can be linked into a user's Geant4
application in order to provide phonon and charge carrier transport in
crystals.  Users must reference G4CMP in their application build in order to
utilize these features.  

After running one of the setup scripts mentioned above (`g4cmp_env.csh` or
`g4cmp_env.sh`), several environment variables will be defined to support
linking G4CMP into your applications:

| Environment variable | Meaning              |  Value in Make build | Value in CMake build |
| ---------------------| ---------------------| ---------------------|----------------------|
| G4CMPINSTALL | Path to g4cmp_env.* scripts  | <path-to-G4CMP> | $CMAKE_INSTALL_PREFIX/share/G4CMP |
| G4CMPLIB | Directory containing libG4cmp.so | $G4WORKDIR/lib/$G4SYSTEM | $G4CMPINSTALL/lib |
| G4CMPINCLUDE | Path to library/include      | $G4INSTALL/library/include | $CMAKE_INSTALL_PREFIX/include |
| G4LATTICEDATA | Path to CrytalMaps directory | $G4INSTALL/CrystalMaps | $G4INSTALL/CrystalMaps |

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


## Application Examples

In addition to the library, G4CMP is distributed with an `examples`
directory containing three simple applications to demonstrate features of
the library.  

* The `phonon` example shows phonon transport and scattering, including
downconversion and mode mixing, in a cylindrical crystal.

* The `charge` example shows electron and hole transport with NTL ("Luke")
emission of phonons and intervalley scattering.

* The `sensor` example shows how to configure the geometry to collect and
record phonon energy by absorption on superconducting TES-style surface
sensors.

Users may copy any of the individual example directories to their own work
area and adapt them as necessary, or use them as inspiration in developing a
more complex experimental model application.

### Building Examples In Situ

If the G4CMP libraries are being built with Make, any of the three
demonstration programs (phonon, charge) may be built as a normal GEANT4 user
application directly from the package top-level directory.  Use the command

	make examples

to build them all, or

	make <name>

to build just one (where <name> is the directory name of interest).  The
executables will be named `g4cmpPhonon` and `g4cmpCharge`, respectively, and
will be written to `$G4WORKDIR/bin/$G4SYSTEM/`.

### Building Examples With CMake

Each example has been set up as a standalone "project" for CMake.  Copy the
example directory, and use `cmake` with `-D` options to set up and build the
example.


## Versioning Information

G4CMP provides somewhat limited access to version information at run time.

Since the package is primarily distributed through GitHub, users can query
the state of their local clone at the command line, using `git describe` to
get back a string such as "G4CMP-190" or "g4cmp-V07-02-02".

For static (tar-ball) distributions, the Git state at the time the tar-ball
was created (using the `make dist` target) will be stored in a file named
`.g4cmp-version`.  This same file will be created as part of the build
process using either Make or CMake (see above).

At runtime, the version string will be available through a call to
`G4CMPConfigManager::Version()`.


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
| **InterValley scattering with matrix elements** |
| epsilon | e/e0      | Relative permittivity     |                    |
| neutDens | N        | Number density of neutron impurities | /volume |
| alpha   |  val      | Non-parabolicity of valleys | energy^-1 (/eV)  |
| acDeform | val      | Acoustic deformation potential | energy (eV)   |
| ivDeform | val val ... | Optical deformation potentials | eV/cm      |
| ivEnergy | val val ... | Optical phonon thresholds     | energy (eV) |
| **InterValley scattering  (Linear and Quadratic Models) ** |
| ivModel     | name | IVRate (matrix), Linear or Quadratic   | string |
| ivLinRate0  | val | Constant term in linear IV expression   | Hz     |
| ivLinRate1  | val | Linear term in linear IV expression     | Hz     |
| ivLinPower  | exp | Exponent: rate = Rate0 + Rate1* E^exp   | none   |
| ivQuadRate  | val | Coefficient for quadratic IV expression | Hz     |
| ivQuadField | val | Minimum field for quadratic IV expression | V/m  |
| ivQuadPower | exp | Exponent: rate = Rate*(E^2-Field^2)^(exp/2) | none |


## Surface Interactions

Transport of both phonons and charge carriers will involve interactions at
the surface of a crystal volume.  The "boundary processes" are modelled on
Geant4's optical physics process, and support reflection, transmission (from
one lattice-equipped volume to another), and absorption with configurable
probabilities.

User applications should use the `G4CMPSurfaceProperty` class, or an
application-specific subclass.  This class has `G4MaterialPropertiesTable`
objects for phonons and charges seaprately; the base class constructor takes
a long list of arguments to fill those tables with common parameters:

  G4CMPSurfaceProperty(const G4String& name,
                       G4double qAbsProb, // Prob. to absorb charge carrier
                       G4double qReflProb, // If not absorbed, prob to reflect
                       G4double eMinK, //Min wave number to absorb electron
                       G4double hMinK, //Min wave number to absorb hole
                       G4double pAbsProb, // Prob. to absorb phonon
                       G4double pReflProb, // If not absorbed, prob to reflect
                       G4double pSpecProb, //Prob. of specular reflection
                       G4double pMinK, //Min wave number to absorb phonon
                       G4SurfaceType stype = dielectric_dielectric);

These parameters are sufficient to model absorption or reflection of both
charges and phonons at the surface of a crystal.  User applications may
choose to define both skin surfaces (for bare crystal substrates) and border
surfaces (with associated sensor/device volumes attached to the crystal)
with different property parameters.

User applications with active sensors for either phonons or charges (or
both), should define a subclass of `G4CMPVElectrodePattern` for each of
those sensors.  If the sensors require additional parameters, those should
be assigned to the material properties table that goes with the surface
above.  See below for a discussion of `G4CMPPhononElectrode`.

Phonon sensors typically involve a superconducting film to couple the
substrate to a sensor (SQUID, TES, etc.).  The `G4CMPKaplanQP` class
provides a parametric model for that coupling, implmenting Kaplan's model
for energy exchange between phonons and quasiparticles from broken Cooper
pairs.  This class expects to find the following material properties defined
for the metal film (defined using the function
`G4MaterialPropertyTable::SetConstProperty("key", value);`).

| Property Key     | Definition                  | Example value (Al) |
|------------------|-----------------------------|--------------------|
| filmThickness    | Thickness of film           | 600.*nm            |
| gapEnergy        | Bandgap of film material    | 173.715e-6*eV      |
| lowQPLimit       | Minimum bandgap multiple for quasiparticles | 3. |
| phononLifetime   | Phonon lifetime in film at 2*bandgap | 242.*ps   |
| phononLifetimeSlope | Lifetime dependence vs. energy | 0.29         |
| vSound           | Speed of sound in film      | 3.26*km/s          |
| subgapAbsorption | Probability to absorb energy below 2*bandgap | 0. |

The last parameter is optional.  It only applies if there is a sensor
involved which is sensitive to heat energy, in which case phonons below
2.*bandgap energy should be treated as directly absorbed with the specified
probability.

A concrete "electrode" class, `G4CMPPhononElectrode`, is provided for simple
access to `G4CMPKaplanQP` from user applications.  An instance of
`G4CMPPhononElectrode` should be registered to the `G4CMPSurfaceProperty`
associated with the phonon sensors' surface.  The material properties listed
above should be registered into the surface's material property table, via
`G4CMPSurfaceProperty::GetPhononMaterialPropertiesTablePointer()`; this
table will be passed into `G4CMPKaplanQP` automatically when it is
registered.  

`G4CMPPhononElectrode` also supports an additional material property,
"filmAbsorption", to specify the "conversion efficiency" for phonons
incident on the registered sensor.  This assumes that the sensor is
implemented as a dedicated volume with an associated border surface.  If
individual sensor shapes are not implemented, this parameter may also
include geometric coverage.
