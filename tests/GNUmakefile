# G4CMP/phonon/GNUmakefile -- for building phonon propagation demo
# $Id: bb6e2ce10a89515fa367d15dcff13f8b03e7a4c3 $
#
# 20160609  Support different executables by looking at target name
# 20170923  Add testChargeCloud
# 20220921  G4CMP-319 -- Add testTemperature
# 20221104  G4CMP-340 -- Move phononKinematics to tools/ directory
# 20250102  G4CMP-436 -- Add testNRyield to exercise Lindhard (NIEL) functions
# 20250428  G4CMP-465 -- Add testSolidUtils for validating transforms in class.

TESTS := electron_Epv latticeVecs luke_dist testBlockData testCrystalGroup \
	g4cmpEFieldTest testChargeCloud testPartition testNRyield \
	testHVtransform testFanoFactor testTemperature testSolidUtils

.PHONY : $(TESTS)

ifndef G4CMP_NAME
help :			# First target, in case user just types "make"
	@echo "G4CMP/tests : This directory contains standalone unit tests"
	@echo
	@echo "luke_dist        : Generate tab-delimited file of NTL phonons"
	@echo "electron_Epv     : Generate tab-delimited file of e- kinematics"
	@echo "latticeVecs      : Apply lattice and valley rotations to vectors"
	@echo "testBlockData    : Demonstrate use of data container"
	@echo "testCrystalGroup : Validate non-orthogonal crystal axes"
	@echo "g4cmpEFieldTest  : Validate COMSOL field file in rectangular box"
	@echo "testChargeCloud  : Validate performance of G4CMPChargeCloud"
	@echo "testHVtransform  : Check lattice transforms and inversions"
	@echo "testFanoFactor   : Verify Fano fluctuations given mean, F"
	@echo "testTemperature  : Exercise thermal distribution functions"
	@echo "testNRyield      : Exercise Lindhard yield (NIEL) functions"
  @echo "testSolidUtils   : Validate the transforms in the SolidUtils class"
	@echo
	@echo Please specify which one to build as your make target, or \"all\"

all : $(TESTS)

$(TESTS) :
	@$(MAKE) G4CMP_NAME=$@ bin

clean :
	@for t in $(TESTS) ; do $(MAKE) G4CMP_NAME=$$t clean; done
else
include $(G4CMPINSTALL)/g4cmp.gmk
endif
