# G4CMP/library/GNUmakefile	Build common-code library libG4cmp.(a|so|dylib)
# $Id$

# Follow GEANT4 library building convention

name := G4cmp

ifndef G4INSTALL
$(error Please configure your GEANT4 environment before building G4CMP)
endif

include $(G4INSTALL)/config/architecture.gmk

# Dependencies of all the phonon classes
CPPFLAGS += -I$(G4BASE)/global/management/include \
	-I$(G4BASE)/global/HEPgeometry/include \
	-I$(G4BASE)/geometry/management/include \
	-I$(G4BASE)/geometry/volumes/include \
	-I$(G4BASE)/particles/management/include \
	-I$(G4BASE)/processes/management/include

# Library goes into user's G4WORKDIR, not into system install
G4LIBDIR := $(G4WORKDIR)/$(G4SYSTEM)

include $(G4INSTALL)/config/common.gmk
