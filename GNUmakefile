# G4CMP/GNUmakefile	Top-level driver to build library and demos
# $Id$

.PHONY : library phonon charge channeling	# Targets named for directory

# Initial target provides guidance if user tries bare |make|
help :
	@echo "G4CMP/GNUmakefile: Drives building library and demos" ;\
	 echo "User must have configured their environment for GEANT4," ;\
	 echo "using the geant4make.sh or .csh script." ;\
	 echo ;\
	 echo "Targets available:" ;\
	 echo "all           Builds everything: library and examples" ;\
	 echo "lib, library  Builds library, including qhull" ;\
	 echo "qhull         Builds qhull general convex hull package" ;\
	 echo "examples      Builds all examples, plus library if needed" ;\
	 echo "phonon        Builds pure phonon example" ;\
	 echo "charge        Builds charge-carrier (e-/h) example" ;\
	 echo "channeling    Builds charged-particle lattice guiding" ;\
	 echo "clean         Remove libraries and examples" ;\
	 echo ;\
	 echo "Users may pass targets through to directories as well:" ;\
	 echo "  make <dir>.<target>"

# User targets

all : lib examples
lib : library
examples : phonon charge channeling

clean : library.clean examples.clean

# Package dependencies

library : qhull
phonon charge channeling : library

# Directory targets

qhull :
	-$(MAKE) -C qhull-2012.1 DESTDIR=$(G4WORKDIR) \
	  BINDIR=$(G4WORKDIR)/bin/$(G4SYSTEM) \
	  LIBDIR=$(G4WORKDIR)/lib/$(G4SYSTEM) \
	  SO=$(if $(findstring Darwin,$(G4SYSTEM)),dylib,$(if $(findstring Win,$(G4SYSTEM)),dll,so))

library phonon charge channeling :
	-$(MAKE) -C $@

library.% \
phonon.% \
charge.% \
channeling.% :
	-$(MAKE) -C $(basename $@) $(subst .,,$(suffix $@))

qhull.% :
	-$(MAKE) -C qhull-2012.1 $(subst .,,$(suffix $@))

# FIXME: This should work with dependencies, but doesn't
examples.% :
	-$(MAKE) phonon.$(subst .,,$(suffix $@))
	-$(MAKE) charge.$(subst .,,$(suffix $@))
	-$(MAKE) channeling.$(subst .,,$(suffix $@))
