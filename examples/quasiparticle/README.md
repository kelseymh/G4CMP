# Quasiparticle Propagation in G4CMP

This example is a demonstration of several features and physics lists that are new to G4CMP as of November 2025. Ryan Linehan, linehan3@fnal.gov made this, so please email him with questions/compliments/complaints.

## Physics Introduction

The physics involved in this example is an extension of the basic physics in G4CMP, and includes several new processes.

* Phonon transmission through surfaces
	*  The physics involved in transmitting phonons through the boundary follows the same logic in `G4CMPBoundaryUtils::ApplyBoundaryAction`, but now in the case that all of absorption-at-electrode (i.e. KaplanQP), "simple" absorption, and reflection fail to trigger, the phonon will transmit through the interface. You can test this by turning off KaplanQP and setting `pAbsProb` and `pReflProb` both to values less than 1.0.
 	*  Boundary surfaces must be applied in both directions for a given interface!
  	*  Currently, no "physics" is done at these interfaces -- phonons pass straight through if they pass through. Proper phonon refraction based on acoustic properties, etc, is an ongoing project.
  	*  Specular vs. diffuse reflection is handled as it was in the previous version.
  	*  Notably, phonons in thin films travel in the full 3-dimensional space (to be contrasted with QPs, next)
*  Phonon Polycrystalline Grain Boundary Scattering: `G4CMPPhononPolycrystalElasticScattering.cc`
	*  This is just an elastic scattering that redirects the phonon after drawing a next step based on a characteristic length. The length represents the characteristic grain boundary size in a polycrystal.
* Bonafide tracking of Bogoliubov Quasiparticles (QPs) in thin-film superconducting volumes
	*  All physics of QPs in thin superconducting films is purely in two dimensions. Right now, the two dimensions are XY, but this will be generalized soon.
* Quasiparticle boundary interactions: `G4CMPQPBoundaryProcess.cc`
	*  `G4CMPSurfaceProperties` now has two additional parameters: `qpAbsProb` (QP Absorption Probability) and `qpReflProb`, which come after the charge and phonon values.
 	*  When you define a superconducting volume (see below) you will have to define a superconducting gap value for that volume's LatticePhysical. If a BogoliubovQP impinges upon a superconductor whose gap is larger than the QP's energy, it will reflect with 100% probability, ignoring the `qpReflProb` you set. If the QP energy is above the gap of the new superconductor, then you transmit with 100% probability unless you have specified a nonzero `qpAbsProb` or `qpReflProb`.
* Cooper-pair breaking by phonons: `G4CMPSCPairBreakingProcess.cc`
  	* The rate of this is given by a dedicated rate function, and is dictated by a combination of parameters set in the `CrystalMaps/Al/config.txt` file and parameters passed into the LatticePhysical attached to a volume. See more in the "Defining Superconductors" section below.
  	* This will produce two BogoliubovQP particles from a phonon above 2*delta
* Phonon radiation by QPs: `G4CMPQPRadiatesPhononProcess.cc`
 	* This will radiate phonons from QPs above delta. The rate is affected by the superconductor parameters disccussed below. 
* QP Recombination: `G4CMPQPRecombinationProcess.cc`
	* This will take a QP and "recombine" it with an ambient quasiparticle that is implicitly in the environment due to some ambient density. A phonon will emerge half of the time, to conserve energy.
 	* This does *not* do n^2 recombination. This recombination is linear in the density of quasiparticles and is a good approximation in the limit of low density of QPs. We'll put back-of-the-envelope numbers to this regime soon. Again, this does *not* do n^2 recombination.
 	* This rate is again dependent on the superconductor parameters discussed below. 
* QP Local Trapping: `G4CMPQPLocalTrappingProcess.cc`
 	* This is a generic linear loss term that kills QPs after they exist for some characteristic lifetime. Notionally this is from trapping on shallow trapping sites, and does not generate phonons.
   	* The rate of this is dependent on a dedicated singular superconductor parameter (see below).
* QP Diffusion: `G4CMPQPDiffusion.cc`
 	* This is a doozy of a function. It uses an efficient MC approach to diffusion in a generalized geometry called Walk-on-Spheres to do diffusion steps of QPs in thin films. Currently only implemented in 2D, and moreover only currently implemented in XY specifically. Will expand to direction agnostic form in a future release.
  	* For fine geometries (like coplanar waveguides), this will take some time to run. The execution time is dependent on the relationship between typical length scales traveled before hitting a boundary and the overall lifetime of the QP (either via recombination, absorption, or local trapping).
  	* If you intend to have QPs in your simulation, this must be turned on for anything to be accurate.

In the G4 macro used to run events, these processes can be turned on and off (mostly) at your discretion. An example code block looks like so:


```
#-------------------------------------
# Turn on/off new processes
#/process/inactivate phononScattering
#/process/inactivate phononDownconversion
/process/inactivate phononPolycrystalElasticScattering
#/process/inactivate qpRecombination
#/process/inactivate qpRadiatesPhonon
#/process/inactivate scPairBreaking
#/process/inactivate qpDiffusion #If QPs can exist, diffusion needs to be active
#/process/inactivate qpLocalTrapping
/process/inactivate qpDiffusionTimeStepper
```

In this example, phononScattering and phononDownconversion are the ``old,'' i.e. already-existing phonon processes, and the rest are part of the TrackedFilmResponse branch. 

> [!CAUTION]
> Most combinations of these processes being on/off work, but QP transport will be UNDEFINED if the qpDiffusion process is turned off. In other words, if you want to do ANY quasiparticle dynamics (say, by looking at "just" phonon pairbreaking), QP diffusion will need to be on. If you want to look at "just" QP phonon radiation, QP diffusion will need to be on as well, etc. The only scenario in which qpDiffusion can be inactive is one in which no tracked QPs are expected to be produced in the simulation.

Final note: for completeness, `qpDiffusionTimeStepper` is a "second" way of doing diffusion, but also requires the qpDiffusion process to be on. It's more of a diagnostic tool and actively slows the code down relative to just using `qpDiffusion`, so most people shouldn't need or want to use it. If you're considering this, may be good to check in with Ryan (linehan3@fnal.gov) to see if it's something you're actually wanting to do.


## Defining Superconductors

On top of the "simple" phonon parameters needed to be defined for any material, there are a set of "superconducting" parameters needed to describe the behavior of QPs and phonon-QP interactions. Below is a table showing these parameters. A few critical notes:
* The first two of these parameters, rigidly named `sc_tau0_qp` and `sc_tau0_ph` are defined in a superconducting material's `config.txt` file. This is done because these parameters are more material-intrinsic than the rest, which may vary depending on the thickness, location, purity, etc., of a physical lattice.
* The rest of these parameters, can in principle be different from superconducting volume to volume on a given chip depending on things like thickness, location, etc., and so we require that you set these for each superconductor volume you define. In particular, for each superconductor volume, you must define a `G4LatticePhysical` object, pass it these parameters, and associate it with that volume. These parameters' names are not rigid, and can be whatever you want as long as you get the argument ordering right for the `G4LatticePhysical` constructor.
* With the table we have provided _example_ values, which are approximately useful values for a boilerplate aluminum film.

| Parameter Name  | Description | Location of Definition | Example Value | Processes Affected |
| ------------- | ------------- | ------------- | ------------- | ------------- |
| `sc_tau0_qp`  | Characteristic QP Lifetime | `CrystalMaps/Al/config.txt` | 438 ns | `G4CMPSCPairBreakingProcess.cc`, `G4CMPQPRecombinationProcess.cc`, `G4CMPQPRadiatesPhononProcess.cc` |
| `sc_tau0_ph` | Characteristic Phonon Lifetime  | `CrystalMaps/Al/config.txt` | 0.242 ns | `G4CMPSCPairBreakingProcess.cc`, `G4CMPQPRecombinationProcess.cc`, `G4CMPQPRadiatesPhononProcess.cc` |
| `polycryElScatMFP` | Characteristic Polycrystalline Grain Boundary Scattering Length | Second argument of `G4LatticePhysical` constructor | 30 nm | `G4CMPPhononPolycrystalElasticScattering.cc` |
| `scDelta0` | Zero-Temperature Superconducting Gap, Δ | Third argument of `G4LatticePhysical` constructor | 180 μeV | `G4CMPSCPairBreakingProcess.cc`, `G4CMPQPRecombinationProcess.cc`, `G4CMPQPRadiatesPhononProcess.cc`, `G4CMPQPBoundaryProcess.cc`, `G4CMPQPDiffusion.cc` |
| `scTeff` | Effective Temperature | Fourth argument of `G4LatticePhysical` constructor | 0.2 K | `G4CMPSCPairBreakingProcess.cc`, `G4CMPQPRecombinationProcess.cc`, `G4CMPQPRadiatesPhononProcess.cc`, `G4CMPQPBoundaryProcess.cc`, `G4CMPQPDiffusion.cc` |
| `scDn` | Normal-state QP Diffusion Constant | Fifth argument of `G4LatticePhysical` constructor | 6 μm^{2} / ns | `G4CMPQPDiffusion.cc` |
| `scTauQPTrap` | Characteristic QP Local Trapping Time | Sixth argument of `G4LatticePhysical` constructor | 1 ms | `G4CMPQPLocalTrappingProcess.cc` | 

We note that _in addition to_ this table, you will have to define interfacial properties and boundaries for substrate-superconductor interfaces, as well as for superconductor-superconductor interfaces (for now). See the physics introduction section above for more info.

> [!NOTE]
> These conventions, including the presence of certain slightly-degenerate processes (like G4CMPPhononPolycrystalScattering) and the method by which the last four parameters are passed into the `G4LatticePhysical`, may morph relatively soon as near-term style fixes are implemented after version 0 of TrackedFilmResponse goes live. If you get errors, the best place to check is the classes and constructors themselves. Next-best is to email linehan3@fnal.gov if you have questions.
										 
> [!TIP]
> The _strength_ of the various processes' dependence on the various parameters varies somewhat widely. For example, if the effective temperature parameter `scTeff` is a small (<10%) fraction of the T_{c} for a given superconductor, then QPs will functionally take *forever* for these QPs to recombine. If your execution is taking forever, this could be a culprit.

## Geometry Description

In this example, the geometry is a superconducting qubit chip with six "candlestick" qubits adjacent to a superconducting feedline. Phonons may transmit from the silicon chip into this superconducting layer (as well as the copper housing meant for thermalization), and may produce quasiparticles that diffuse in 2D around the layer.

Coming Soon: Geometry tips/alerts and the way to circumvent the issues

## Primary Event

## Execution and Output

* For now, the step output accessed using the commands in QuasiparticleSteppingAction is only meaningful while executing with a single worker thread -- multithreading is not implemented at the moment.

## General Tips/Notes

* For now, need one G4LatticeLogical and G4LatticePhysical for every physical superconducting volume. This may incur a slight calculational cost at the beginning of a given run, since for each there is a calculation of the QP scattering/recombination taus as a function of energy relative to the gap. The more dedicated volumes you have, the longer you should expect the startup to take.


## Tutorial Walkthrough (Coming soon!)






