# G4CMP Documentation

## Overview

G4CMP is a framework designed to add to the Geant4 toolkit for use in condensed
matter and low-energy physics. Developed for the low-temperature community,
G4CMP is capable of modeling several physics processes relevant to phonon and
charge collection at cryogenic temperatures
[.](../README.md) These include anisotropic phonon
transport and focusing, phonon isotope scattering, anharmonic downconversion,
oblique charge carrier propagation with inter-valley scattering, and emission of
Neganov-Trofimov-Luke phonons by accelerated carriers
[.](https://arxiv.org/pdf/1403.4984.pdf)

[//]: # (For second source: sections 2 and 3)

The package provides a collection of particle types, physics processes, and
supporting utilities for this purpose
[.](../README.md)
It has been used by collaborators at the
SuperCDMS project to successfully reproduce theoretical predictions and
experimental observations such as phonon caustics, heat pulse propagation times,
and mean carrier drift velocities[;](https://arxiv.org/pdf/1403.4984.pdf)
the G4CMP package, however, is sufficiently
general that it is useful for other experiments employing cryogenic phonon and/
or ionization detectors[.](https://arxiv.org/pdf/1403.4984.pdf)

[//]: # (For second source: section 4. For last source: section 1)

## [Release Notes](https://github.com/kelseymh/G4CMP/blob/master/ChangeHistory)

## [Installation](https://github.com/kelseymh/G4CMP#building-the-package)

## [Examples](https://github.com/kelseymh/G4CMP#application-examples)

## Settings and Usage

### [User Environment](https://github.com/kelseymh/G4CMP#user-environment)

### [Defining Crystal Dynamics](https://github.com/kelseymh/G4CMP#defining-the-crystal-dynamics)

### [Surface Interactions](https://github.com/kelseymh/G4CMP#surface-interactions)

## Testing

Not Currently Available

## Contributing

No documentation Currently Available

## Crystal Physics

Some things G4CMP models:

- Acoustic phonons, electrons, and holes in cryogenic crystals.
- Anisotropic phonon propagation, oblique carrier propagation and phonon emission
by accelerated carriers.

### Crystal Lattices

G4CMP comes with configuration data for germanium and silicon crystals in the
package itself using the G4LogicalLattice and G4PhysicalLattice classes; these
are found in the CrystalMaps directory, where users can also define crystals of
other materials.
These config.txt files are plain text with names, values, and units. User
applications must specify the config name separately from the G4Material name.
Each lattice definition requires several sections:

- Crystal parameters
- Phonon parameters
- Charge carrier parameters
- Hole and electron masses

The material properties and crystal structure are implemented via the
G4LogicalLattice class, which provides the natural coordinate frame of the
lattice and associates it to a specific Geant4 “placement volume” with an
orientation. The G4PhysicalLattice class handles local/lattice/valley coordinate
transforms[.](https://wordpress-new.physics.tamu.edu/tobackgroup/wp-content/uploads/sites/21/2022/08/HEP-Software-Foundation-Colloquium-Kelsey.pdf)

[//]: # (For source: slides 20 and 36)

### Charge Transport

Physically, charge transport occurs when incident particles promote electrons to
the conduction band, simultaneously creating holes. However, there are several
peculiarities of charge transport that G4CMP simulates which are worthy of note.

#### Valleys and Intervalley Scattering

The lowest energy bands in crystals have particular orientations, called valleys.
Electrons travel along valleys but are also scattered between them.

Electrons are transported along valleys because they have different effective
masses parallel and perpendicular to the valley axis—this means that they might
have different masses in multiple directions. This is modeled by the Electron
Mass Tensor. In particular, letting the valley axis be $\vec{x}$,

$$
m_= =
\begin{bmatrix}
  m_{\parallel} &0 &0 \\
  0 &m_{\perp} &0 \\
  0 &0 &m_{\perp}
\end{bmatrix}
, \ \vec{p} = m_= \vec{v}, \ E = \vec{p}^T
m_=^-1 \vec{p}
$$

This relationship only applies close to the valley axis, and in general the Mass
tensor is direction-dependent. G4CMP uses a fixed mass tensor for all kinematics.

Intervalley scattering occurs when electrons are strongly scattered by
absorption of thermal phonons—a large momentum transfer can move an electron
from one valley to another. This scattering contributes to electron drift
speed.
G4CMP models this process with the G4CMPInterValleyScattering class, with
several options for determining the scattering rate. For that purpose, the use
may use G4CMPIVRateLiner, G4CMPIVRateQuadratic, or G4CMPInterValleyRate.

[//]: # (ibid, slides 21, 37-38)

#### Phonon Emission

Charges accelerated in an electric field can radiate phonons. These Phonons are
known as Neganov-Trofimov-Luke (NTL) Phonons, and are emitted when charges
accelerated by an E-field attain velocities well above vsound. This occurs when
the charges interact with the lattice, radiating phonons, reducing their energy,
and changing direction. The phonon emission rate can be modeled by,

$$
\nu = \frac{3l_0}{v_{sound}} \frac{Ma}{(Ma-1)^3}
$$

With total phonon emission equalling energy gained from the potential difference
across which the charge accelerates. Here, $l_0$ is the scattering length for
the charge type and $Ma$ is the Mach number ($v/v_{sound}$) for the charge.

The class G4CMPLukeScattering implements this phenomenon, with rate determined
by the G4CMPLukeEmissionRate class.

[//]: # (ibid, slide 39)

#### Charge Recombination

One “fact of life” in the Geant4 framework is that particles are independent and
isolated—they do not mutually interact. However, in the low-energy case, charges
(i.e. electrons and hotels) at the surfaces of a crystal cannot escape because
they do not “reflect” against the bias voltage. As a solution, G4CMP assumes
that these charges recombine with some pre-existing partner (an electron with a
hole or a hole with an electron).

When this occurs, half of the bandgap energy is released as phonons at the
material’s Debye frequency.

- Silicon: 15THz, 62.03 meV
- Germanium: 2THz, 8.27 meV

Since electron-hole pairs are created initially, both charges recombining and
releasing half the bandgap energy ensures energy conservation.

This process is modeled in G4CMP in the G4CMPDriftRecombinationProcess class.

[//]: # (ibid, slide 40)

#### Charge Trapping on Impurities

Similarly to recombination, this occurs when charges become stuck in place.
Trapping, however, occurs when charges are stopped by impurities in the bulk of
the detector material. When particles are trapped at a shallow (~ meV) depth,
the bandgap energy is not recovered.

There are two types of impurities in the crystal lattice, resulting in four
kinds of capture, one for each kind of charge and impurity:

- $ e + D^0 → D^—, \ e + A^+ → A^0, \ h + A^0 → A^+, \ h + D^— → D^0 $

These captures have separate rates for electrons and holes. They can also be
device-dependent and even history (neutralization) dependent.

Charges trapped in this way can contribute to the detector’s charge collection
signal if they are trapped near enough to the electrodes.

This process is modeled with the G4CMPDriftTrappingProcess class, with the
electron and hole trapping lengths managed by the corresponding configuration
commands:
    /g4cmp/electronTrappingLength
    /g4cmp/holeTrappingLength
These parameters are global for an entire computing job; that is, G4CMP assumes
that the user is simulating a single device. Versions of these parameters that
are attachable to individual detector volumes may be implemented in the future.

[//]: # (ibid, slide 41)

#### Impurity Trap Reionization

Trap Reionization is essentially the inverse of trapping. Tracks can interact
with traps that have immobilized a charge, releasing charges. When this occurs
at a shallow (~meV) depth, the bandgap energy is not absorbed.

Two types of impurities here result in four types of reionization each with a
separate rate.

- $ e + D^− → 2 e + D^0 , \ e + A^+ → e + h + A^0 $
- $ h + A^0 → 2 h + A^+ , \ h + D^− → h + e + D^0 $

As before, rates can depend both on device and neutralization history.

The G4CMPDriftTrapIonization class handles this process in the package, and has
corresponding configuration commands.
    /g4cmp/eDTrapIonizationMFP
    /g4cmp/eATrapIonizationMFP
    /g4cmp/hDTrapIonizationMFP
    /g4cmp/hATrapIonizationMFP
As in the previous section, these are global parameters.

[//]: # (ibid, slide 42)

### Phonons

Phonons are a type of energy-carrying quasiparticle; in particular, they consist
of quantized lattice oscillations that occur in several ways. Phonons can either
be longitudinal (compression waves), or transverse (shear waves), and can
propagate in either low energy (“acoustic”) or high energy (“optical”) states.

[//]: # (ibid, slide 22)

#### Phonon Mode Group Velocity

Using the crystal stiffness matrix along a given direction $\hat{n}$, we have
the Christoffel matrix
$ D_{il} = C_{ijlm} \cdot \hat{n}^j \cdot \hat{n}^m / \rho $,
whose eigenmodes are phase velocity and polarization. Group velocity is
then computed from these factors.

In order to speed up processing, G4CMP generates lookup tables with steps of
$ \hat{n} $ coordinates and interpolates between steps. These processes are
governed by the G4CMPPhononKinematics class.

[//]: # (ibid, slide 44)

#### Phonon Impurity Scattering and Anharmonic Decay

Phonons can scatter off of impurities in the crystal lattice, changing their
mode, from longitudinal to slow or fast transverse, for example. The rate of
this scattering scales like E4. Specifically, $ \nu = B \cdot (E/h)^4 $,
where $ B = 2.43 \cdot 10^{-42} s^3 $ in Silicon.

G4CMP implements this phenomenon using wavevector (energy) conservation. A
different mode is chosen based on the configured density of states, and uses the
corresponding wavevector to determine the phonon’s new velocity vector. This is
governed by G4PhononScattering and supplemented by G4CMPPhononScatteringRate.

[//]: # (ibid, slides 23 and 45)

Another possibility after scattering, besides transforming modes, is splitting
into pairs of various modes. Longitudinal (L) phonons can do this, splitting
either into two transverse (T) phonons or a new longitudinal (L’) and transverse
(T) phonon. The rate of this process scales like $ E^5 $. In particular,
$ \nu = D \cdot (E/h)^5 $, where $ D = 2.43 \cdot 10^{-42}s^3 $ and the fraction
of decays to TT compared to L’T pairs is 0.74 in Silicon.

The splitting process equipartitions early “hot” (Debye energy in the tens of
meV) phonons into a sea of meV-scale phonons. G4CMP uses G4PhononDownConversion
and G4CMPDowncoversionRate classes to manage the simulation of this process.

After early high-energy phonons split in the wake of an energy deposit, the
detector crystal becomes filled with a “gas” of low-energy (≲ meV) phonons with
all modes represented, moving in all directions. Sensors on the top and bottom
of the crystal can absorb phonons to measure the magnitude of the energy deposit.

[//]: # (ibid, slides 23, 46-47)

### Energy Partitioning in G4CMP

Geant4 typically doesn’t produce “trackable” electrons below tens of eV.
Instead, it records an “energy deposit” value associated with the electron’s
parent track. In this method, $ dE/dx $ summarizes all the conduction electrons
produced by a track. However, in semiconducting crystals, the minimum energy
required to generate one electron-hole pair is the bandgap, at around 1 eV. In
the end, the typical pair energy is 3-4 eV per pair, with some variation.

Further, ions (including alpha particles) induce motion in nearby atoms in the
lattice. This results in Non-ionizing energy loss (NIEL) and athermal phonons,
each with a Debye energy in the tens of meV.

G4CMP addresses these issues via the G4CMPSecondaryProduction and
G4CMPEnergyPartition classes, allowing for charge tracking at the lower energy
levels used in low-temperature measurement that is otherwise impossible in
Geant4.

[//]: # (ibid, slide 48)

The Relative magnitude of $ dE/dx $ versus NIEL for ions depends on the charge
and mass of the projectile as well as atomic number and mass of the crystal
atoms. The ionization yield is computed by taking $ dE/dx $ as a fraction of the
total. G4CMP does this in its code for ion hits, but it is now done
automatically in Geant4 10.7: the issue is addressed via forward-compatibility,
because G4CMP will not recalculate the yield in the case of non-zero NIEL. This
process remains in the G4CMPLindhardNIEL and G4CMPLewinSmithNIEL classes.

[//]: # (ibid, slide 49)

## QET Physics

G4CMP is primarily used for simulations in the material of a detector crystal,
but it also supports some methods and classes for simulating detector response
to events within the crystal. Here, we cover some of its applications in that
realm.

### Kaplan Quasiparticle Model

This model is used for thin-film superconductors common in cryogenic electronic
sensors, such as TESs or QETs, KIDs, and qubits. These sensors are populated
with Cooper pairs, which can be broken upon phonon absorption into electron
quasiparticles. For small films, quasiparticles transport faster than the
thermal response. As for Crystal properties, aluminum film parameters are
supplied with the package via a Geant4 properties table, including film
thickness, Cooper pair gap energy (2Δ), $ v_{sound} $, and phonon lifetime.

G4CMP models energy transfer, QP transport, and phonon re-emission via the
G4CMPKaplanQP class.

[//]: # (ibid, slide 50)

### G4CMPKaplanQP

In G4CMPKaplanQP, the simulation of the Kaplan Quasiparticle Model is treated as
instantaneous; it iterates to find the equilibrium state and contains no
time-dependent information.

Substrate phonon absorption in the thin film is modeled by using the mean free
path and thickness to find the probability of this occurring. In particular,

$$
P=e^{\frac{-4d}{MFP}}, \ MFP = \frac{v_{sound}}{\tau(E)}, \
\tau(E)=\frac{\tau_0}{1+\delta \tau(E/\Delta-2)}
$$

where P is the probability of substrate phonon absorption, d is the film
thickness, MFP is the mean free path, and 𝛕 is the phonon lifetime.

Phonon energy in these simulations goes to break Cooper pairs, and becomes the
QP energy. QP or phonon energy is absorbed into the tungsten TES, at which time
some QP energy goes back into phonons via QP “decay” (emission), and some phonon
energy is re-emitted back into the substrate. In G4CMPKaplanQP, the processing
loop ends when the available phonon energy is zero. This model is not suitable
for “bare” films that lack an attached energy absorber.

[//]: # (ibid, slide 51)

### Phonon Readout Model

Phonon energy deposits are collected in time bins resulting in a matching
readout. Coupled differential equations model the electrothermal response of
TESes, the bias current, indicative (SQUID) coupling, and other phenomena. CVODE
(from LLNL) is used to solve for current output in each time bin, while
configuration files specify detector components and characteristics including
heat flow, resistances, inductances, TESes per channel, and others.

[//]: # (ibid, slide 52)
