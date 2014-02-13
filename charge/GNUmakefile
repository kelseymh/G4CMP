# G4CMP/charge/GNUmakefile -- for building charge-carrier propagation demo
# $Id$

name := g4cmpCharge
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
$(error Please configure your GEANT4 environment before building G4CMP)
endif

.PHONY: all
all: lib bin

EXTRALIBS := -L$(G4WORKDIR)/lib/$(G4SYSTEM) -lG4cmp -lqhullcpp

include $(G4INSTALL)/config/binmake.gmk
CPPFLAGS += -I../library/include -I$(G4WORKDIR)/include -Dqh_QHpointer

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*

