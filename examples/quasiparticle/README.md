

 =================================================================
                  Quasiparticle propagation in G4CMP
 =================================================================
                         Ryan Linehan - FNAL
                          linehan3@fnal.gov

This example demonstrates an integration of phonon transmission through interfaces as well, phonon-quasiparticle interactions, and quasiparticle dynamics as well.

1.INTRODUCTION

The physics involved in this example is an extension of the basic physics in G4CMP, and includes several new processes.

* Phonon transmission through interfaces
  	* Now, `G4CMPSurfaceProperties` must be defined with an additional two parameters: `qpAbsProb` (QP Absorption probability) and `qpReflProb`, but once these are defined, you can enable passage of phonons through a surface by setting `pAbsProb` and `pReflProb` to values less than 1.
  	* Boundary surfaces must be applied in both directions for a given interface!
	* Currently, no "physics" is done at these interfaces -- phonons pass straight through if they pass through. Specular vs. diffuse reflection is handled as it was in the previous version.
* Cooper-pair breaking by phonons: `G4CMPSCPairBreakingProcess.cc`
  	* The rate of this is given by a dedicated rate function, but is ultimately dictated by parameters set in the `CrystalMaps/Al/config.txt` file. Only a few of them should be tweaked by a user at the moment: one is an effective temperature `sc_Teff`, which should be less than Tc for a given crystal.
  	* This will produce two BogoliubovQP particles from a phonon above 2*delta
* Phonon radiation by QPs: `G4CMPBogoliubovQPRadiatesPhononProcess.cc`
 	* This will radiate phonons from QPs above delta. The rate is again dictated by that effective temperature in `CrystalMaps/Al/config.txt`.
* QP Recombination: `G4CMPBogoliubovQPRecombinationProcess.cc`
	* This will take a QP and "recombine" it with an ambient quasiparticle that is implicitly in the environment due to some ambient density. A phonon will emerge half of the time, to conserve energy.
 	* This does *not* do n^2 recombination. This recombination is linear in the density of quasiparticles and is a good approximation in the limit of low density of QPs. We'll put back-of-the-envelope numbers to this regime soon. Again, this does *not* do n^2 recombination.
 	* This rate is *strongly* dependent on the Teff you use in `CrystalMaps/Al/config.txt`. If you set this to below about 10% of Tc for a given superconductor, you be waiting *forever* for these QPs to recombine. 
* QP Local Trapping: `G4CMPBogoliubovQPLocalTrappingProcess.cc`
 	* This is a generic linear loss term that kills QPs after they exist for some characteristic lifetime. Notionally this is from trapping on shallow trapping sites
   	* This is another crystal parameter, `sc_tau_qptrap` at the moment, located in `CrystalMaps/Al/config.txt`.
* QP Diffusion: `G4CMPBogoliubovRandomWalkTransport.cc`
 	* This is a doozy of a function. It uses an efficient MC approach to diffusion in a generalized geometry called Walk-on-Spheres to do diffusion steps of QPs in thin films. Currently only implemented in 2D, and moreover only currently implemented in XY specifically. Will expand to direction agnostic form in a future release.
  	* For fine geometries (like coplanar waveguides), this will take some time to run. The execution time is dependent on the relationship between typical length scales traveled before hitting a boundary and the overall lifetime of the QP (either via recombination, absorption, or local trapping).
  	* If you intend to have QPs in your simulation, this must be turned on for anything to be accurate.
* Gap Engineering: `G4CMPBogoliubovQPRandomWalkBoundary.cc`
 	* QPs can also move between superconducting volumes that are all in-plane, but are prevented from entering a superconductor whose gap is higher than the QP's energy.

This example will simulate an application in which all of these physics processes are playing an active role in the evolution of the system.

2. GEOMETRY

In this example, the geometry is a superconducting qubit chip with four "candlestick" qubits adjacent to a superconducting feedline. Phonons may transmit from the silicon chip into this superconducting layer (as well as the copper housing meant for thermalization), and may produce quasiparticles that diffuse in 2D around the layer.

3. PHYSICS PROCESSES

In the G4 macro used to run events, there are a set of new processes that you can turn on and off during your simulations. An example code block looks like so:

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

In this example, phononScattering and phononDownconversion are the ``old,'' i.e. already-existing phonon processes, and the rest are part of the TrackedFilmResponse branch. Most combinations of these processes being on/off work, but a critical thing to keep in mind is that QP transport will be UNDEFINED if the qpDiffusion process is turned off. In other words, if you want to do ANY quasiparticle dynamics (say, by looking at "just" phonon pairbreaking), QP diffusion will need to be on. If you want to look at "just" QP phonon radiation, QP diffusion will need to be on, etc. The only scenario in which qpDiffusion can be inactive is one in which no QPs are expected to be produced in the simulation.

Final note: qpDiffusionTimeStepper is a "second" way of doing diffusion, but also requires the qpDiffusion process to be on. It's more of a diagnostic tool and actively slows the code down relative to just using qpDiffusion, so most people shouldn't need or want to use it. If you're considering this, may be good to check in with Ryan (linehan3@fnal.gov) to see if it's something you're actually wanting to do.

4. PRIMARY EVENT

TBD

5. EXECUTION & OUTPUT

TBD

6. TESTING

TBD

7. GENERAL TIPS

* For now, need one G4LatticeLogical and G4LatticePhysical for every physical superconducting volume. This may incur a slight calculational cost at the beginning of a given run, since for each there is a calculation of the QP scattering/recombination taus as a function of energy relative to the gap. The more dedicated volumes you have, the longer you should expect the startup to take.
* For now, need to define boundaries between not only substrate and thin film, but also between different volumes in the thin film. This includes boundaries between mother and daughter volumes
* Right now, the CrystalMaps config.txt files are only effective:
 	* For Al, the "boilerplate" phonon parameters are equivalent to that of Si. We'll update these to accurate values for Al in the next week or so. The superconducting parameters (the last 8 parameters in the file) are the ones that contain the new physics.
	* For Nb, the parameters are basically equivalent to Al, which is factually incorrect. For the moment, the difference is effectively just Teff, Tc, and the energy gap. Would not really recommend using Nb for now, until we flesh this out better.
  



