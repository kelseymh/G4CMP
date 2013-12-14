# G4CMP/GNUmakefile	Top-level driver to build library and demos
# $Id$

.PHONY : library phonon charge channeling	# Targets named for directory

# Initial target provides guidance if user tries bare |make|
help :
	@echo "G4CMP/GNUmakefile: Drives building library and demos"
	@echo "User must have configured their environment for GEANT4,"
	@echo "using the geant4make.sh or .csh script."
	@echo
	@echo "Targets available:"
	@echo "all     	Builds everything: library and examples"
	@echo "lib     	Builds library"
	@echo "examples	Builds examples only, library only if necessary"
	@echo "phonon  	Builds pure phonon example"
	@echo "charge  	Builds charge-carrier (e-/h) example"
	@echo "channeling"
	@echo "        	Builds example of charged-particle lattice guide"
	@echo "clean   	Remove libraries and examples"
	@echo
	@echo "Users may pass targets through to directories as well:"
	@echo "  make <dir>.<target>"

# User targets

all : lib examples
lib : library
examples : phonon charge channeling

library phonon charge channeling :
	-$(MAKE) -C $@

clean : library.clean examples.clean

# Directory targets

library.% \
phonon.% \
charge.% \
channeling.% :
	-$(MAKE) -C $(basename $@) $(subst .,,$(suffix $@))

# FIXME: This should work with dependencies, but doesn't
examples.% :
	-$(MAKE) phonon.$(subst .,,$(suffix $@))
	-$(MAKE) charge.$(subst .,,$(suffix $@))
	-$(MAKE) channeling.$(subst .,,$(suffix $@))
