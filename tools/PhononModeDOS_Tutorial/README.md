# Tutorial for Calculation of Phonon DOS and Band Structure

For computing the phonon DOS and band structure necessary for implementing materials in G4CMP, preexisting DFT engines and python packages
can be used. Quantum Espresso (QE) offers Open-Source codes optimized for performing electronic-structure and nanoscale materials modeling. 
The code utilizes density-functional theory, psuedopotentials and plane waves to perform the computation. Phonopy is a python package that 
utilizes DFT engines, like QE, to compute phonon related properties for a given crystal.

##Step 0: Configure the environment

Install QE: https://www.quantum-espresso.org/Doc/user_guide/node7.html
Install Phonopy: https://phonopy.github.io/phonopy/install.html

##Step 1: Define your Crystal

This example will follow calculation of the phonon DOS and band structure for beta Sn. 

Once you know the material and crystal properties, find the crystal of interest in the Materials Project material database. Download the 
CIF and open the file to validate the final structure is correct. 

For beta Sn, the material can be found here: https://legacy.materialsproject.org/materials/mp-84/# 

Downloading the Conventional Standard CIF, one can confirm that the lattice parameter are appropriate for their configuration (a=b=5.91,
c=3.25)

##Step 2: Prepare the QE input file

In order to begin our computation, first we must perform a self-consistant field calculation to find the electron density and effective 
potential. This is an iterative process that guesses the initial density, builds an effective potential, solves the Kohn-Sham equations, 
computes the electron density, then compares to computed density to the initial density until the values agree within a defined 
threshold. To perform this calculation, we must create our input file using the CIF.

Materials Project provides a generator for the input file based on the CIF file: https://tools.materialscloud.org/qeinputgenerator/ 

Once this is generated, you should have a file like this as well as corresponding pseudopotentials:


	&CONTROL
  	  calculation = 'scf'				
  	  etot_conv_thr =   4.0000000000d-05  ->  Convergence threshold on total energy (a.u) for ionic minimization
   	  forc_conv_thr =   1.0000000000d-04  ->  Convergence threshold on forces (a.u) for ionic minimization
  	  outdir = './out/' ->  Directory where output files are stored 
  	  prefix = 'Sn'  ->  Prefix for generated files
  	  pseudo_dir = './pseudo/'  ->  Directory where the psuedopotential is stored  
  	  tprnfor = .true.  ->  calculate forces
  	  tstress = .true.  ->  calculate stress
  	  verbosity = 'high'
	/
	&SYSTEM
  	  degauss =   1.4699723600d-02  ->  value of the gaussian spreading (Ry) for brillouin-zone integration in metals
  	  ecutrho =   4.8000000000d+02  ->  Kinetic energy cutoff (Ry) for charge density and potential
  	  ecutwfc =   6.0000000000d+01  ->  kinetic energy cutoff (Ry) for wavefunctions
  	  ibrav = 0  ->  Bravais-lattice index. Optional only if space_group is set.
  	  nat = 4  ->  Number of atoms
  	  nosym = .false.  ->  if set to true, symmetry is not used
  	  ntyp = 1  ->  Number of different elements
  	  occupations = 'smearing'  ->  gaussian smearing for metals
  	  smearing = 'cold'  ->  Marzari-Vanderbilt-DeVita-Payne cold smearing (see PRL 82, 3296 (1999))
	/
	&ELECTRONS
  	  conv_thr =   8.0000000000d-10  ->  Convergence threshold for SCF calculation
  	  mixing_beta =   4.0000000000d-01  ->  mixing factor for self-consistency
	/

	CELL_PARAMETERS {angstrom}  ->  Defines a, b and c
	5.9080080000 0.0000000000 0.0000000000
	0.0000000000 5.9080080000 0.0000000000
	0.0000000000 0.0000000000 3.2456220000

	ATOMIC_SPECIES -> Atom  Atomic Weight  Pseudopotential File
	Sn 118.71 sn_pbesol_v1.4.uspp.F.UPF

	ATOMIC_POSITIONS {crystal}  ->  Atomic positions matching the CIF
	Sn 0.0000000000 0.0000000000 0.0000000000 
	Sn 0.5000000000 0.0000000000 0.2500000000 
	Sn 0.5000000000 0.5000000000 0.5000000000 
	Sn 0.0000000000 0.5000000000 0.7500000000 


Take care in setting the convergence threshold. It should be d-10 or lower to ensure an accurate calculation.

##Step 3: Run the SCF calculation

Now that we have the input file, we can begin the SCF calculation. First we must create our ideal and displaced supercells:

	% phonopy --qe -d --dim="2 2 2" -c pwscf.in

where '2 2 2' corresponds to the number of displaced supercells created. For this command, you should have three new files. Supercell-001.in
and supercell-002.in correspond to displaced supercells, while supercell.in describes a perfect supercell.

To allow us to use the executable from QE, we must add a header specifying parameters needed by QE:

 	&control
       	   calculation = 'scf'
    	   tprnfor = .true.
    	   tstress = .true.
    	   pseudo_dir = './pseudo/'
    	   disk_io = 'none'
 	/
 	&system
    	   ibrav = 0
    	   nat = 32  -> ensure this corresponds to the number of atoms in your supercell files
    	   ntyp = 1
    	   ecutwfc = 70.0
 	/
 	&electrons
    	   diagonalization = 'david'
     	   conv_thr = 1.0d-10
 	/
	K_POINTS automatic
 	2 2 2  1 1 1

This is saved in the header.in file and can be added with the following command:

	% for i in {001,002};do cat header.in supercell-$i.in >| betaSn-$i.in; done

Now we perform the SCF calculation for both displaced cells:

	% pw.x -i betaSn-001.in |& tee betaSn-001.out
	% pw.x -i betaSn-002.in |& tee betaSn-002.out

Each command took 1hr to run on a single core.

##Step 4: Calculate atomic forces

Using our results from Step 3, we can then compute the FORCE_SETS with the following command:

	% phonopy -f betaSn-001.out betaSn-002.out

##Step 5: Calculate the phonon band structure

To perform this calculation, we must create a short input file:

	DIM = 2 2 2
	BAND = 0.0 0.0 0.0  0.5 0.0 0.0  0.5 0.5 0.0  0.0 0.0 0.0  0.5 0.5 0.5

Naming this file band.conf, we can then compute the band structure with the following command:

	% phonopy --qe -c pwscf.in -p band.conf

##Step 6: Prepare files for extraction of mode dependent DOS

The file we need to plot the mode dependent DOS is mesh.yaml. This file is genrated from phonopy_disp.yaml and phonopy_params.yaml.
These can be created from our previous computations. 

	% phonopy --qe -d --dim 2 2 2 --pa auto -c pwscf.in  ->  phonopy_disp.yaml
	% phonopy --sp -f betaSn-00{1,2}.out  ->  phonopy_params.yaml

Now we can create mesh.yaml:

	% phonopy-load --band "0.0 0.0 0.0  0.5 0.0 0.0  0.5 0.5 0.0  0.0 0.0 0.0  0.5 0.5 0.5" --mesh 41 41 41 --pdos "1, 2" -p

##Step 7: Post-processing for G4CMP

G4CMP_Extract.py extracts and plots the mode dependent DOS contributions for the first three acoustic modes. The relative DOS should
be extracted at 1 THz for use with G4CMP. Additionally, the Debye energy input needed should correspond to the maximum acoustic phonon
energy (maximum of mode 3).
