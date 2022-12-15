# G4CMP/GNUmakefile	Top-level driver to build library and demos
# $Id: b50f1e17e36e73b8ee4a2e6a8e10327f3864cbbd $
#
# Add Mac and Windows handling for Qhull build (we aren't using CMake there)
# Add "dist" target to make tar-ball for distribution to non-SLAC users
# Temporarily exclude "channeling" from all examples, until it builds
# Add G4CMP_SET_ELECTRON_MASS support to build with dynamic mass code
# Move individual examples to "examples/" directory, add tools and tests
# Move QHull down to library/ directory; update tar-ball generator
# Drop G4CMP_SET_ELECTRON_MASS code blocks; not physical
# Add new "sensors" example directory
# Add top-level ".g4cmp-version" file to track as-built tag/SHA
# Add Geant4 version checking
# Manually set version with G4CMP_VERSION=xxx if Git not available
# Add pass-through of thread-safety "code sanitizer" flags
# Split XXX.% targets to ensure everything gets built properly

# G4CMP requires Geant4 10.4 or later
g4min := 10.4

.PHONY : library phonon charge tests tools	# Targets named for directory
.PHONY : all lib dist clean qhull examples

# Initial target provides guidance if user tries bare |make|
help :
	@echo "G4CMP/GNUmakefile: Drives building library and demos" ;\
	 echo "User must have configured their environment for GEANT4" ;\
	 echo "$(g4min) or later, using the geant4make.sh or .csh script." ;\
	 echo ;\
	 echo "Targets available:" ;\
	 echo "all           Builds everything: library and examples" ;\
	 echo "lib, library  Builds library, including qhull" ;\
	 echo "examples      Builds all examples, plus library if needed" ;\
	 echo "phonon        Builds pure phonon example" ;\
	 echo "charge        Builds charge-carrier (e-/h) example" ;\
	 echo "sensors       Builds FET digitization sensor example" ;\
	 echo "tools         Builds support utilities (lookup table maker)" ;\
	 echo "tests         Builds small test programs for classes" ;\
	 echo "clean         Remove libraries and examples" ;\
	 echo ;\
	 echo "Users may pass targets through to directories as well:" ;\
	 echo "  make <dir>.<target>" ;\
	 echo ;\
	 echo "For developers ONLY, make a distribution tar-ball:" ;\
	 echo "dist          Builds a tar ball of the code, excluding .git/" ;\
	 echo ;\
	 echo "For step-by-step debugging, set G4CMP_DEBUG=1";\
	 echo "For thread-safety checking, set G4CMP_USE_SANITIZER=1"

# User targets

all : lib examples tests tools
lib : library
examples : phonon charge sensors

clean :		# FIXME: This doesn't work as dependencies
	-$(MAKE) tests.clean
	-$(MAKE) tools.clean
	-$(MAKE) examples.clean
	-$(MAKE) library.clean 

dist : g4cmp.tgz

# Version identification file (found under .../share in CMake build

G4CMP_VERSION_FILE := .g4cmp-version

.PHONY : version
version :
	@([ -d .git ] && git describe || echo $(G4CMP_VERSION)) \
	   > $(G4CMP_VERSION_FILE)

# Directory targets

export G4CMP_DEBUG		# Turns on debugging output
export G4CMP_USE_SANITIZER	# Enables multithread checking
export G4CMP_SANITIZER_TYPE

library :
	-$(MAKE) version
	-$(MAKE) -C $@

library.% :
	-$(MAKE) -C $(basename $@) $(subst .,,$(suffix $@))

tests.% :
	-$(MAKE) -C $(basename $@) $(subst .,,$(suffix $@))

tools.% :
	-$(MAKE) -C $(basename $@) $(subst .,,$(suffix $@))

phonon charge sensors : library
	-@$(MAKE) -C examples/$@

phonon.% \
charge.% \
sensors.% :
	-$(MAKE) -C examples/$(basename $@) $(subst .,,$(suffix $@))

# FIXME: These should work with dependencies, but don't
examples.% :
	-$(MAKE) phonon.$(subst .,,$(suffix $@))
	-$(MAKE) charge.$(subst .,,$(suffix $@))
	-$(MAKE) sensors.$(subst .,,$(suffix $@))

tests : tests.all
tools : tools.all

# Make source code distribution (construct using symlinks and tar -h)

g4cmp.tgz : clean
	-$(MAKE) version
	@/bin/rm -rf G4CMP ;\
	 mkdir G4CMP ;\
	 ln -s ../README.md ../LICENSE G4CMP ;\
	 ln -s ../CMakeLists.txt ../cmake G4CMP ;\
	 ln -s ../GNUmakefile ../g4cmp.gmk G4CMP ;\
	 ln -s ../g4cmp_env.sh ../g4cmp_env.csh G4CMP ;\
	 ln -s ../G4CMPOrdParamTable.txt G4CMP ;\
	 ln -s ../library ../examples ../tests ../tools G4CMP ;\
	 ln -s ../CrystalMaps G4CMP ;\
	 ln -s  ../$(G4CMP_VERSION) G4CMP ;\
	 gtar -hzc -f $@ G4CMP ;\
	 /bin/rm -rf G4CMP

# Check for minimum Geant4 version (set at top of file)

g4ver := $(shell geant4-config --version)
ifneq ($(g4min), $(firstword $(sort $(g4min) $(g4ver))))
  $(error Geant4 $(g4min) or later required.  Using $(g4ver))
endif
