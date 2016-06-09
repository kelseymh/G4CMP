# G4CMP/channeling/GNUmakefile -- for building channeling propagation demo
# $Id$

name := g4cmpChannel
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
$(error Please configure your GEANT4 environment before building G4CMP)
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*

