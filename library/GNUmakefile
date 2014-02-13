# G4CMP/library/GNUmakefile	Build common-code library libG4cmp.(a|so|dylib)
# $Id$

# Follow GEANT4 library building convention

name := G4cmp

ifndef G4INSTALL
$(error Please configure your GEANT4 environment before building G4CMP)
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -I$(G4INCLUDE)	# Always defined in 9.6 and later

# Library goes into user's G4WORKDIR, not into system install
G4LIBDIR := $(G4WORKDIR)/lib/$(G4SYSTEM)

include $(G4INSTALL)/config/common.gmk
