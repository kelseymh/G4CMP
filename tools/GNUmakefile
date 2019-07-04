# G4CMP/phonon/GNUmakefile -- for building phonon propagation demo
# $Id: bb6e2ce10a89515fa367d15dcff13f8b03e7a4c3 $
#
# 20160518  Use G4CMPINSTALL instead of ".." to find includes
# 20160609  Support different executables by looking at target name

# Add additional utility programs to list below
TOOLS := g4cmpKVtables
.PHONY : $(TOOLS)

ifndef G4CMP_NAME
help :			# First target, in case user just types "make"
	@echo "G4CMP/tools : This directory contains standalone utilities"
	@echo
	@echo "g4cmpKVtables : Generate phonon K-Vgroup mapping files"
	@echo
	@echo Please specify which one to build as your make target, or "all"

all : $(TOOLS)

$(TOOLS) :
	@$(MAKE) G4CMP_NAME=$@ bin

clean :
	@for t in $(TOOLS) ; do $(MAKE) G4CMP_NAME=$$t clean; done
else
include $(G4CMPINSTALL)/g4cmp.gmk
endif
