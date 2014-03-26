# G4CMP/GNUmakefile	Top-level driver to build library and demos
# $Id$
#
# Add Mac and Windows handling for QHull build (we aren't using CMake there)
# Add "dist" target to make tar-ball for distribution to non-SLAC users
# Temporarily exclude "channeling" from all examples, until it builds

.PHONY : library phonon charge channeling	# Targets named for directory
.PHONY : all lib dist clean qhull examples

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
	 echo "  make <dir>.<target>" ;\
	 echo ;\
	 echo "For developers ONLY, make a distribution tar-ball:" ;\
	 echo "dist          Builds a tar ball of the code, excluding .git/" ;\
	 echo ;\
	 echo "For step-by-step debugging, set G4CMP_DEBUG=1"

# User targets

all : lib examples
lib : library
examples : phonon charge ### channeling

clean : library.clean examples.clean qhull.cleanall

dist : g4cmp.tgz

# Package dependencies

library : qhull
phonon charge channeling : library

# Directory targets

export G4CMP_DEBUG	# Turns on debugging output

library phonon charge channeling :
	-$(MAKE) -C $@

library.% \
phonon.% \
charge.% \
channeling.% :
	-$(MAKE) -C $(basename $@) $(subst .,,$(suffix $@))

# Manually configure building the QHull libraries in Geant4 style
ISMAC := $(findstring Darwin,$(G4SYSTEM))
ISWIN := $(findstring Win,$(G4SYSTEM))
DYLIB_OPTS := -dynamiclib -undefined suppress -flat_namespace

qhull :
	-$(MAKE) -C qhull-2012.1 DESTDIR=$(G4WORKDIR) \
	  BINDIR=$(G4WORKDIR)/bin/$(G4SYSTEM) \
	  LIBDIR=$(G4WORKDIR)/lib/$(G4SYSTEM) \
	  CC_OPTS3="$(if $(ISMAC),$(DYLIB_OPTS),)" \
	  SO=$(if $(ISMAC),dylib,$(if $(ISWIN),dll,so)).6.3.1 \
	  all install

qhull.% :
	-$(MAKE) -C qhull-2012.1 $(subst .,,$(suffix $@))

# FIXME: This should work with dependencies, but doesn't
examples.% :
	-$(MAKE) phonon.$(subst .,,$(suffix $@))
	-$(MAKE) charge.$(subst .,,$(suffix $@))
	### -$(MAKE) channeling.$(subst .,,$(suffix $@))

# Make source code distribution (construct using symlinks and tar -h)

g4cmp.tgz : clean
	@mkdir G4CMP ;\
	 ln -s ../README ../GNUmakefile ../g4cmp.gmk G4CMP ;\
	 ln -s ../library ../qhull-2012.1 ../charge ../phonon G4CMP ;\
	 ln -s ../channeling ../CrystalMaps G4CMP ;\
	 gtar -hzc --exclude '.*' -f $@ G4CMP ;\
	 /bin/rm -rf G4CMP
