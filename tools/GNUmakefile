# G4CMP/phonon/GNUmakefile -- for building phonon propagation demo
# $Id: bb6e2ce10a89515fa367d15dcff13f8b03e7a4c3 $
#
# 20160518  Use G4CMPINSTALL instead of ".." to find includes
# 20160609  Support different executables by looking at target name
# 20221104  G4CMP-340 -- Move phononKinematics and plotting utility here.

# Add additional utility programs to list below
TOOLS := g4cmpKVtables phononKinematics
.PHONY : $(TOOLS) plot_phonon_kinematics.py


ifndef G4CMP_NAME
help :			# First target, in case user just types "make"
	@echo "G4CMP/tools : This directory contains standalone utilities"
	@echo
	@echo "g4cmpKVtables : Generate phonon K-Vgroup mapping files"
	@echo "phononKinematics : Generate phonon kinematics and plot"
	@echo
	@echo Please specify which one to build as your make target, or \"all\"

all : $(TOOLS) plot_phonon_kinematics.py

$(TOOLS) :
	@$(MAKE) G4CMP_NAME=$@ bin

plot_phonon_kinematics.py : 
	@/bin/cp -f $@ $(G4WORKDIR)/bin/$(G4SYSTEM)/$@

clean :
	@for t in $(TOOLS) ; do $(MAKE) G4CMP_NAME=$$t clean; done
	@/bin/rm f $(G4WORKDIR)/bin/$(G4SYSTEM)/plot_phonon_kinematics.py
else
include $(G4CMPINSTALL)/g4cmp.gmk
endif
