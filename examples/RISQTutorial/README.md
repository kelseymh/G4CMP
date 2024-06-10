# RISQ G4CMP Advanced Tutorial 

Assembled by Ryan Linehan, linehan3@fnal.gov

This tutorial is meant to provide a somewhat lengthy introduction to how to use G4CMP for a variety of applications. We focus here on a geometry and a set of analyses that are relevant to the superconducting qubits (QIS) field. 

## Preliminaries

### Installing Geant4 and G4CMP
While we leave a thorough discussion of the installation procedures to the main G4CMP readme file, it is useful to have a short reminder of this to establish some directory names that we'll use throughout the rest of the tutorial. We'll start with a reminder that in order to run this example, you'll need to install ROOT and both the geant4 and G4CMP packages. On my machine, each of these has three directories associated with its build: a source directory `XXXXX`, a build directory `XXXXX-build`, and an install directory `XXXXX-install`. On my machine, the base name (`XXXXX`) for my geant4 build is `geant4.10.07.p04`, and the base name for the G4CMP build is `G4CMP_RISQTutorial`. NB: the last steps in the process of each of these installations should be to run `make` and `make install` while in the `XXXXX-build` directory, for both geant4 and G4CMP


> [!IMPORTANT]
> To streamline your ability to prep for this tutorial, we recommend installing Geant4 with the cmake flags `-DGEANT4_INSTALL_DATA=ON` and `-DGEANT4_USE_OPENGL_X11=ON`, and to build with C++14. In particular, the OpenGL flag will enable visualization, which we will frequently use. However, if you can successfully run other visualizers like DAWN, those are also perfectly fine.

### Setting up environment
Assuming you've built these directories and you're opening up a new terminal, you'll need to source the environmental setup scripts for these:
```
source /path/to/geant4.10.07.p04-install/bin/geant4.sh
source /path/to/G4CMP-install/share/G4CMP/g4cmp_env.sh
```
Now we can make our example. Copy this tutorial's source directory into a new directory -- I like to copy it outside of the whole G4CMP source directory just to avoid confusion and remember that this is its own executable that needs to be made. Moreover, make build and install directories to accompany it:
```
cd /path/to/G4CMP_RISQTutorial 
cp -r ./examples/RISQTutorial /path/to/
cd /path/to/
mkdir RISQTutorial-build
mkdir RISQTutorial-install
```
Now we head into our build directory and run CMake:
```
cd RISQTutorial-build
cmake -DCMAKE_INSTALL_PREFIX=/path/to/RISQTutorial-install -DCMAKE_CXX_STANDARD=14 ../RISQTutorial/
```
If this runs successfully, we should be able to run make and then make install, and we're done:
```
make
make install
```
If those build without errors, we should be ready to get started.



## Example 1: Superconducting Qubit Chip Geometry
Let's begin our tutorial by understanding how to build a somewhat complicated detector geometry. We'll put together a chiop geometry based on a specific type of transmon superconducting qubit design called an Xmon (paper [here](https://arxiv.org/pdf/1304.2322)), which has a cross-shaped island coupled to a quarter-wave resonator. Our chip geometry is derived from designs made by the McDermott group at UW Madison, and has the following features:
* A silicon chip onto which a superconducting ground plane is laid
* Six transmon qubits with accompanying resonators, patterned into the ground plane.
* Six control lines, one per qubit, with a wirebond pad, patterned into the ground plane
* A central transmission line, with in/output wirebond pads, patterned into the ground plane
* Chip housing, through which the chip is thermalized.
To see these features, let's run our first macro: `vis.mac`. Since we are already in our `/path/to/RISQTutorial-build` directory, we can just run the executable:
```
./RISQTutorial
```
This should boot up G4CMP into an interactive session, from which you can execute either single commands or, as we do here, macros. Run the following to run `vis.mac`.
```
/control/execute ../RISQTutorial/G4Macros/vis.mac
```
What you will ideally see (contingent upon which viewing software you have installed) is something like the following:
<img width="400" alt="FirstVis" src="https://github.com/kelseymh/G4CMP/assets/20506221/4b3cc3e1-7a7d-418c-ba86-a2877d414914">


Here we can see the various features of the qubit chip. The chip outline itself is shown in gray, the superconducting features are shown in cyan, and the thermal mount (copper) is shown in orange.

In this example, we discuss some strategies for building relatively complex geometries like this one. We start this discussion within the `src/RISQTutorialDetectorConstruction.cc` file, which is the hub for geometry-building in this example. Within this, let's take a look within the `SetupGeometry()` function. While we will not go into gory detail within this function (since that discussion can be found in other examples), a few important points should be made:
* The silicon chip, copper housing, and ground plane are all established with the world volume as their logical mother volume.
* The transmission line, resonator assemblies (including qubits), and flux lines are all constructed using dedicated classes. These classes can be found in the `src` and `include` directories. Constructing these objects in classes enables tiling of these complex geometries across a chip as we do for the resonator assemblies. 
* The transmission line, resonator assemblies (including qubits), and flux lines are all established with the ground plane as their logical mother volume.

Navigating to the `src/RISQTutorialResonatorAssembly.cc` file, we see that the resonator assembly is made of three elements: a "base layer," a "resonator line," and a "shunt capacitor cross" (i.e. the qubit). The base layer is a block layer into which the other two elements are placed. The other two elements are constructed from central superconducting volumes nested inside vacuum volumes, to form, for example, the coplanar waveguide making the resonator. Some other important notes:
* It is important to give all volumes and sub-volumes unique names -- this allows us to be able to later access phonon hits in specific volumes when running phonon simulations. If you decide to use multiple copies of the custom objects like resonator assemblies, it is therefore important to pass in a unique identifier to each. Here, we use six copies of the ResonatorAssembly object but pass in a unique name for each copy in the `pName` parameter.
* Some classes, like the `src/RISQTutorialTransmissionLine.cc` are made from other sub-classes like `src/RISQTutorialPad.cc`.
* Each class has a member object `fFundamentalVolumeList`, which is a list of fundemental/irreducible sub-volumes within the class object. This volume list is used to iteratively establish boundary surfaces between the silicon substrate, which exists only in the top-level `src/RISQTutorialDetectorConstruction.cc` file, and the volumes created within the dedicated class object. The calls to these fundamental sub-volumes are done in the `src/RISQTutorialDetectorConstruction.cc` file. The surface properties themselves are also defined in the `SetupGeometry()` function within `src/RISQTutorialDetectorConstruction.cc`. We will discuss these in the following example.

Other assemblies, such as `src/RISQTutorialTransmissionLine.cc` follow a similar organizational strategy as the resonator assembly, with the caveat that we have constructed the envelopes (spatial boundaries) of all assemblies in such a way that no assembly envelopes overlap. As the ground plane is being used as the mother volume for all of these assemblies, ensuring this lack of overlaps is useful for uniquely defining a volume for every point on the superconductor film.

It is useful to have one place to define most or all hardcoded detector geometry constants. We use a file `include/RISQTutorialDetectorParameters.hh`, which we include in all `.cc` files that define geometry. This detector parameters file defines dimensions for all sub-volumes of every detector class, but the relative size and positioning (defined in terms of those dimensions) of volumes is done within the class files.


## Example 2: Simple Phonon Collection Efficiency Study
Now that we've fleshed out the geometry construction a bit, let's do a study of phonon collection efficiency in this chip. 

### Throwing a phonon
Let's start out first by just simulating a single phonon and understanding the dynamics. Navigate to your `RISQTutorial-build` directory and execute the following command:
```
./RISQTutorial
/control/execute ../RISQTutorial/G4Macros/throwPhonon.mac
```
This should yield a visualization that looks a bit like this:

<img width="400" alt="ThrowPhonon_pabs10pct" src="https://github.com/kelseymh/G4CMP/assets/20506221/80273712-9a95-4191-aa47-2ee3e9d63160">


What happened in this event? We threw a single longitudinal acoustic phonon of 30 meV from a point midway between the bottom and top of the chip (here at 4.62 mm and 5 mm in Z, respectively), and watched it undergo downconversion, isotopic scattering, and surface scattering. For more detailed information, let's take a look at one of the tracking snippets printed out: 

<img width="953" alt="ThrowPhononFirstPrintout" src="https://github.com/kelseymh/G4CMP/assets/20506221/136de436-5a0e-422a-9e82-e535ba9fd9fe">


The first meaningful step (Step 1) of our initial longitudinal phonon (phononL) is downconversion, which produces a fast transverse phonon (phononTF) and a slow transverse phonon (phononTS). In the visualization, phononL/phononTF/phononTS are shown in blue/green/red, respectively. A snippet of the phononTS track is also visible: while it cannot itself undergo downconversion, it can undergo isotopic scattering, which it does several times. The transportation step visible is one in which this phonon reflects off the bottom of the chip. We also take a look at the phononTF track produced:

<img width="959" alt="ThrowPhononLaterInPrintout" src="https://github.com/kelseymh/G4CMP/assets/20506221/c0a3ed76-ba4e-43ce-813c-318e821584fc">


Here, we notice something similar to the phononTS, but with a downconversion step at the end. Transverse phonons are not allowed to downconvert in G4CMP, the isotopic scattering process enables mode-mixing between phonon polarizations, which is exactly what is happening here: our phononTF changes modes via scattering until it reaches a phononL which then decays prior to another mode change. While this does change the polarization of the phonon under the hood (making the physics correct), this is not easily reflected in the tracking printout.)

Now let's change some of the parameters in our chip. One of the most important parameters to get comfortable with is the phonon absorption/reflection coefficients at interfaces. We modify this by opening up the `src/RISQTutorialDetectorConstruction.cc` file and looking for the line
```
fSiNbInterface = new G4CMPSurfaceProperty("SiNbInterface",
                                          1.0,0.0,0.0,0.0,
                                          0.1,1.0,0.0,0.0);
```
This line governs much of the physics at the boundary between the silicon substrate and the superconductor. The first four entries after the interface name deal with e/h absorption and reflection -- let's ignore those for now. The last four entries deal with phonon absorption and reflection. In order, the first three entries of those four indicate that we have a 10% probability of phonon absorption at the surface, a 100% probability of reflection if that absorption fails, and if reflection happens, it will happen specularly with a 0% chance.

> [!IMPORTANT]
> It's generally a good idea to read through the base G4CMP code to understand how the code works, but it is particularly useful to read through the code governing how boundary parameters you set actually lead to reflection and absorption. The logic blocks governing actions of quanta at boundaries is in the `G4CMPBoundaryUtils::ApplyBoundaryAction()` function as shown below:
> <img width="629" alt="image" src="https://github.com/kelseymh/G4CMP/assets/20506221/1c6b4ae5-5b14-4ad2-910a-92f280679f73">


Let's change our absorption probability from 0.1 to 0.01, to allow our phonon to ricochet around the detector for longer. After doing this, we need to re-make and then rerun:
```
cd /path/to/RISQTutorial-build
make
make install
./RISQTutorial
/control/execute ../RISQTutorial/G4Macros/throwPhonon.mac
```
which should yield something like the following visualization:

<img width="400" alt="ThrowPhonon_pabs1pct" src="https://github.com/kelseymh/G4CMP/assets/20506221/9abd1a7d-a9a3-48c2-9189-4bedc6022021">


Phonons now propagate farther in our chip, which makes sense: every time they hit the ground plane or any superconducting structures (which covers much of the top side of the chip) they now only have a 1% of absorption. 

What's the right value to use for an absorption probability? This bit of microphysics is tricky and benefits from experimentally probing a given chip, but a reasonable starting point for a first simulation is an analytical expression [here](https://arxiv.org/abs/2404.04423) (with variable descriptions at the link):

<img width="550" alt="image" src="https://github.com/relineha/TestRISQTutorial/assets/20506221/1b66c33f-6783-443b-a8df-69011f356157">

In short, a naive absorption probability depends on the thickness of the superconducting film, its gap, the phonon energy, and a few other parameters specific to the film material.

Before proceeding, let's do two more things. First, let's take a look at the output of the simulation. This output for this tutorial is contained within two files: `RISQTutorial_primary.txt` and `RISQTutorial_hits.txt`. The former should currently only contain a single line showing information about the single primary phonon thrown: 
<img width="930" alt="image" src="https://github.com/kelseymh/G4CMP/assets/20506221/2ca178a3-1e9c-4ab4-8812-742b54d7ef39">

The latter should contain the set of "hits" that occur, where phonons are absorbed: 

<img width="930" alt="image" src="https://github.com/kelseymh/G4CMP/assets/20506221/789b0751-d5e0-4c78-b6e4-cdb07e3d8129">


In this event we had six phonons create hits on the superconductor and thermal bath mounts on the underside of the chip. One useful thing to note is that the energies deposited in these hits are at the few-to-several meV scale. This reflects the characteristic energies of phonons in the ballistic regime in silicon, where boundary scatters occur with significantly higher rates than further downconversion.

> [!TIP]
> Homework problem: Saving this output in a set of text files isn't the most space-efficient way to store the data we're generating, and it's kind of annoying to have two output files for the same set of events. Can you think of a better way to do this?

The second thing we should do before proceeding to the phonon collection efficiency study is to go back into our detector construction file and set our phonon absorption probability to 0.5. As you might discover, phonons that are allowed to bounce around longer will result in a longer runtime for a given simulation. For the purpose of getting reasonably high stats for this tutorial, we want to limit that a bit. Go ahead and change that back go 0.5, and then remake the example:
```
cd /path/to/RISQTutorial-build
make
make install
```

### Phonon Collection Efficiency Study
Now let's actually do some chip characterization. We'll use a new macro for this: `pceStudy.mac`. Some notable differences between this macro and `throwPhonon.mac` are:
* `/vis/` commands (removed): Since we're going to throw lots of phonons, all of the visualization and trajectory information is removed.
* `/gps/energy 0.004 eV`: We're simulating our phonons with different energies: 4 meV instead of 30 meV. This starts the phonons in the ballistic regime, and while downconversion can still happen, these phonons should largely not downconvert in this simulation.
* `/gps/pos/` commands: We're simulating our phonons with a different spatial configuration, with bou. Here the events are simulated uniformly throughout the chip volume rather than at a single point. This is useful so that we can get a sense of how the phonon collection efficiency varies as a function of position in the chip.
* `/run/beamOn 200000`: We're now simulating 200k single-phonon events, rather than a single phonon. 
* `/tracking/verbose 0`: Since we're running so many events, we're going to turn the tracking verbosity off.
* `/g4cmp/phononBounces 1000`: We're going to increase the number of phonon reflections that are allowed before killing a phonon -- this is 100 by default but we'll extend it a bit.


Let's go ahead and run this using the following commands:
```
cd ./path/to/RISQTutorial-build
./RISQTutorial ../RISQTutorial/G4Macros/pceStudy.mac
```
If we now look in the output files `RISQTutorial_primary.txt` and `RISQTutorial_hits.txt` we see that there are 200k entries in the former and approximately 200k hits in the latter.

Let's start to analyze these events. We have a dedicated ROOT macro that is designed to do this in the `RISQTutorial/AnalysisTools` directory. Execute the following commands:
```
cd /path/to/RISQTutorial/AnalysisTools/
root -l
```
Once ROOT opens up, execute
```
.L RISQTutorialAnalysis.cc
PCEStudy("/path/to/RISQTutorial-build/RISQTutorial_primary.txt","/path/to/RISQTutorial-build/RISQTutorial_hits.txt")
```
When this code is done, it should spit out a file called `AnalysisOutput.root`. Let's take a look inside this with our ROOT TBrowser to look at the histograms filled. First, we have a set of cross-checks on the locations of our primary phonons, so that we can ensure we generated the starting distribution correctly:

<img width="1067" alt="PrimariesCrossCheck200k" src="https://github.com/kelseymh/G4CMP/assets/20506221/1955ed1f-6b40-480f-8ee1-67900387c27f">


These look like they were spawned properly in our chip: they're uniform in X, Y, and Z. Now let's do a cross-check on where our hits are found.

<img width="1178" alt="HitsCrossCheck200k" src="https://github.com/kelseymh/G4CMP/assets/20506221/560a93f7-ba56-4f53-86a5-3c84e4994acc">


Again, things make sense: phonons are absorbed and create hits across the ground plane as well as on the bottom of the chip where it sits on the copper mount to make contact with the thermal bath. We don't see hits in a few locations, specifically where we have gaps between our wirebond pads and the ground plane. In reality, there are similar regions where our CPWs and qubits are located, but they are too small of features to see given our 200k statistics and this histogram binning.

For a phonon collection efficiency study, we don't actually care too much about phonons that hit superconductor far from the qubit islands -- any QPs out in the ground plane will just diffuse around until recombination, producing little to no effect on the qubits themselves. So we want to redefine how we think about a "hit." We do this in `/path/to/RISQTutorial/src/RISQTutorialSensitivity.cc`. This function is where we define what counts as a "hit." Within the `IsHit()` function, let's change our `selectTargetVolumes` variable to true. This will only consider a phonon absorption a `hit` if it lands on a surface bordering a volume with the string `shuntConductor` in its name -- these are the qubit crosses we've built in our geometry. Go ahead and remake:
```
cd /path/to/RISQTutorial-build
make
make install
```
This time, let's increase the number of phonons we throw to two million by editing the `pceStudy.mac` file with the following line
```
/run/beamOn 2000000
```
Then let's rerun the macro:
```
cd /path/to/RISQTutorial-build
./RISQTutorial ../RISQTutorial/G4Macros/pceStudy.mac
```
...and rerun the root analysis...
```
cd /path/to/RISQTutorial/AnalysisTools/
root -l
.L RISQTutorialAnalysis.cc
PCEStudy("/path/to/RISQTutorial-build/RISQTutorial_primary.txt","/path/to/RISQTutorial-build/RISQTutorial_hits.txt")
```
Now, if we look at our hits histogram, we see that the only hits that are registered are those that land within the six qubit islands present in the geometry:

<img width="400" alt="image" src="https://github.com/kelseymh/G4CMP/assets/20506221/c8534396-019e-4518-97a8-0bc2e6f778be">


With this definition of hits, we can now compute a meaningful phonon collection efficiency for phonons impinging upon the qubit islands. We define efficiency as the fraction of phonon energy emitted that is absorbed into any qubit. To do this as a function of spatial location, we need two histograms: one with total energy emitted as a function of primary XY and one with the total energy absorbed in any qubit as a function of primary XY:

<img width="1181" alt="ComponentsToPCE" src="https://github.com/kelseymh/G4CMP/assets/20506221/7ac6ef58-dba0-4c98-a3f7-f75cdacdd4e4">


From this we can compute a phonon collection efficiency:

<img width="400" alt="image" src="https://github.com/kelseymh/G4CMP/assets/20506221/d3e6658f-fe93-4ed2-81d1-4b6345019e4c">


This makes sense: the phonon collection efficiency for relatively large absorption (50%) at superconducting interfaces will imply that only phonons that are relatively close in XY to a qubit island will actually be able to register an absorption in that island. Overall, the PCE looks to be at the O(0.001) level.

> [!TIP]
> Homework problem: We've made an overall phonon collection efficiency map in this example. Can you make a plot showing the phonon collection efficiency in a single qubit as a function of the phonon spawn point's radial distance from the qubit?

## Example 3: Muon Event
Now that we've got a feel for phonon propagation, let's run an example where we look at the full physics chain of a muon passing through our chip. In these events, the muon will deposit O(several hundred keV) along its path through the chip, producing lots of ionization (electron/hole pairs) and phonons. This example will enable us to explore how G4CMP handles the production and propagation of this ionization, as well as how it produces phonons. We'll use the macro `/path/to/RISQTutorial/G4Macros/throwMuon.mac` to simulate a 4 GeV muon passing through the chip at an angle 14 degrees from the plane of the chip.

> [!IMPORTANT]
> In order to simulate conventional "high-energy" particles like muons in G4CMP where we commonly think about meV-eV scale phenomena, we need to remember to include in our example a physics list that recognizes those. We do this "under the hood" in this example within the RISQTutorial.cc file, and set our main HEP physics list to be FTFP_BERT. See this file for more info.

### Visualization to Build Intuition
First, let's just try to visualize this event to help build our intuition. Let's blindly run the macro without modification in an effort to simulate every G4CMPDriftElectron, G4CMPDriftHole, and phonon, and see what happens:
```
cd /path/to/RISQTutorial-build
./RISTutorial
/control/execute ../path/to/RISQTutorial/G4Macros/throwMuon.mac
```
After a while, it's pretty obvious what's going on here. Even if it weren't for the event printouts slowing things down, simulating all quanta produced in this interaction is computationally expensive. If we may expect a 400 keV energy deposition to ultimately end up as a set of 4 meV ballistic phonons, that's 100 million tracks, each with up to 1000 steps. Because that's absurdly impractical for a single event, G4CMP has a set of tools to downsample quanta generation in the event, reducing the computational cost while maintaining reasonable accuracy in simulation response. The knobs that enable us to do this are in the `throwMuon.mac` macro:
```
/g4cmp/producePhonons 0.0
/g4cmp/sampleLuke 0.0
/g4cmp/produceCharges 0.01
```
Here, the arguments give the fraction of recombination-generated phonons, luke-emission-generated phonons, and charge carriers that get created and tracked. We'll start out with just looking for a percent of the charges. Running the macro again, we find:

<img width="400" alt="BasicMuonOnlyCharges1pct" src="https://github.com/kelseymh/G4CMP/assets/20506221/c02aa606-c4c1-4ef2-9446-98c41fa7245b"> <img width="400" alt="BasicMuonOrthographic" src="https://github.com/kelseymh/G4CMP/assets/20506221/56181435-dfd2-49d8-be03-f149105d5d4c">


Here, the white track is the muon, the pink tracks are the G4CMPDriftElectrons, and the orange tracks are the G4CMPDriftHoles. As the muon is generated in the top right region of the chip, this is the origin point for the ionization. Let's zoom in and look at what's going on near the muon track:

<img width="400" alt="MuonEventWithZoom" src="https://github.com/kelseymh/G4CMP/assets/20506221/3df04bcd-7be0-464c-afbd-64d73d372fe1">


Here the G4CMPDriftElectrons are undergoing occasional direction changes. If we turn on tracking (`/tracking/verbose 1`) and rerun this, we can see why: 

<img width="954" alt="LukeScattering" src="https://github.com/kelseymh/G4CMP/assets/20506221/59fa9b36-00ab-4d18-86d4-7932891ac5e1">


Here, the G4CMPDriftElectrons are undergoing Luke phonon emission and radiating phonons, which occurs when G4CMPDriftElectrons have sufficient energy to travel faster than the speed of sound in the substrate. Since we don't see phonons radiated from these scattering points because of our downsampling, let's play around with the downsampling to see if we can confirm this intuition. We choose these somewhat artificial values to see this Luke emission:

```
/g4cmp/producePhonons 0.0
/g4cmp/sampleLuke 0.1
/g4cmp/produceCharges 0.00001
```
and, zooming in, we can find what we're looking for -- phonon tracks sprouting out of kinks in G4CMPDriftElectron tracks:

<img width="400" alt="LukeScatteringPhonons" src="https://github.com/kelseymh/G4CMP/assets/20506221/94e9414e-0055-47be-a668-14c4189d793a">

Okay, so we've confirmed our intuition: when we enable Luke phonon generation, we see them radiated from G4CMPDriftElectrons and G4CMPDriftHoles. 

> [!TIP]
> Homework problem: Can you figure out which piece of code in the G4CMP package handles the modeling of this Luke emission for high-energy G4CMPDriftElectrons?


Zooming out a bit, we can look back a few visualizations and see that G4CMPDriftElectrons and G4CMPDriftHoles ultimately free stream to the edges of the chip, which happens once they drop below the energy threshold for Luke emission. In G4CMP, their interactions with the surfaces of the chip are handled in a similar way to phonons: there are coarse absorption and reflection parameters that determine any further propagation. For this example, we're going to simplify and set the absorption parameter for these charge carriers to be 100% at these boundaries. In our simulations, these electrons and holes will therefore stop, shed their kinetic energy as phonons at these interfaces, and then recombine, each giving up one half of the bandgap energy as phonons. We therefore expect to see some phonon generation around the edges of the chip as well. We can check this, by turning on the last of our downsampling parameters:
```
/g4cmp/producePhonons 0.01
/g4cmp/sampleLuke 0.01
/g4cmp/produceCharges 0.001
```
which give us something like:

<img width="400" alt="RepresentativeEvent" src="https://github.com/kelseymh/G4CMP/assets/20506221/9964ef06-1b53-4660-b1d3-a84b67249eb7">

Now that we have all of our phonons visually represented, we have arrived at an event display that's a little bit more representative of the phonon response of our chip to a muon event (keeping in mind, of course, that we're still suppressing charge and phonons by massive suppression factors).

### Doing a simple single-event analysis
Let's now see what information we can extract from even a single simulated muon event. Let's make our goal to estimate the total phonon energy absorbed by each of the six qubits on this chip. Let's modify our `throwMuon.mac` macro to use the following downsampling parameters:
```
/g4cmp/producePhonons 1.0
/g4cmp/sampleLuke 1.0
/g4cmp/produceCharges 0.01
```
Since this event dumps about 400 keV into the chip and creates O(100,000) e/h pairs, downsampling to one percent of those still gives us O(1000) e/h pairs, which certainly gives us a good spread of G4CMPDriftElectron and G4CMPDriftHole ending positions (i.e. boundary recombination phonon generation locations) while still limiting the amount of time needed to simulate this event. On my laptop, this event took approximately 45-60 seconds to simulate. At the end, we again get a set of primaries and hits in our output files -- looking in `RISQTutorial_hits.txt` we find that we got something like 524 hits across all six qubits. That's not a whole lot of hits and will keep statistical uncertainty in our final in-qubit energies kind of high, but let's call it fine for now. Let's head back to our analysis tools and fire up that analysis macro:
```
cd /path/to/RISQTutorial/AnalysisTools/
root -l
.L RISQTutorialAnalysis.cc
```
This time in ROOT, we use a different dedicated function (hey, that's convenient!) to do analysis of our muon event. To help this analysis, we need to label our qubits as shown in the figure below to help with understanding our analysis results. We execute our function:
```
AnalyzeMuonEvent("../../RISQTutorial-build/RISQTutorial_primary.txt","../../RISQTutorial-build/RISQTutorial_hits.txt",0.01)
```
where the final argument in the `AnalyzeMuonEvent()` function is the downsampling factor we've applied to the `/g4cmp/produceCharges` command -- it corrects the total in-qubit phonon to more accurately represent what would be present without the downsampling. If we now look in our output root file `AnalysisOutput.root`, we find our desired quantity:

<img width="400" alt="QubitLabels" src="https://github.com/kelseymh/G4CMP/assets/20506221/d2b1960d-9fd0-4b01-89d3-2070f53b7de2"> <img width="400" alt="StandardMuonEventOutput" src="https://github.com/kelseymh/G4CMP/assets/20506221/246f9362-4f14-49fe-bdae-bdf70f44f113">

We see that qubit 2 has the largest in-qubit energy by over an order of magnitude -- this makes sense, as this qubit is almost right over the muon track. If we then change our muon location by tweaking the spawn location in our in our `throwMuon.mac` file:
```
/gps/pos/centre 0.0 0.5 0.56 cm
```
we see these relative in-qubit energies shift accordingly:

<img width="400" alt="CenteredMuonEventQubitLabels" src="https://github.com/kelseymh/G4CMP/assets/20506221/341463fe-00bb-4324-82e0-89dbacd2163e"><img width="400" alt="CenteredMuonEventOutput" src="https://github.com/kelseymh/G4CMP/assets/20506221/09f8f7b0-431a-4692-99a6-d85dab2c97fc">


Here, the muon is ACTUALLY right under a qubit, which leads to a significantly higher excess in Qubit 1's energy over that of other qubits. This illustrates that according to G4CMP, we might be able to back out some coarse position information from a muon event if we can reconstruct the in-qubit energy somehow. This is something we can in fact do! For more information, see [this paper](https://arxiv.org/abs/2404.04423).

> [!TIP]
> Homework problem: We're taking a naive sum of all hit energies over time in the above plots. Can you make a plot of the time series of the hits seen by each qubit? Challenge question: how would you go about scaling an "energy-deposited-vs-time" plot to account for the downsampling?

