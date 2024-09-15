# G4CMP: Geant4 Condensed Matter Physics Library

## Overview

G4CMP is a framework designed to add to the Geant4 toolkit for use in
condensed matter and low-energy physics. Developed for the low-temperature
community, G4CMP is capable of modeling several physics processes relevant
to phonon and charge collection at cryogenic temperatures.  These include
anisotropic phonon transport and focusing, phonon isotope scattering,
anharmonic downconversion, oblique charge carrier propagation with
inter-valley scattering, and emission of Neganov-Trofimov-Luke phonons by
accelerated carriers.  

The published descriptions of G4CMP include a 2014 preprint
([arXiv:1403.4984](https://arxiv.org/pdf/1403.4984.pdf)), and a 2023
publication in NIM A
([doi:10.1016/j.nima.2023.168473](https://doi.org/10.1016/j.nima.2023.168473), [arXiv:2302.05998](https://arxiv.org/abs/2302.05998)).

The package provides a collection of particle types, physics processes, and
supporting utilities for this purpose.  It has been used by collaborators at
the SuperCDMS project to successfully reproduce theoretical predictions and
experimental observations such as phonon caustics, heat pulse propagation
times, and mean carrier drift velocities.  The G4CMP package, however, is
sufficiently general that it is useful for other experiments employing
cryogenic phonon and/ or ionization detectors.
