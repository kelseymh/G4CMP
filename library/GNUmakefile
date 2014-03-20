# G4CMP/library/GNUmakefile	Build common-code library libG4cmp.(a|so|dylib)
# $Id$

# Follow GEANT4 library building convention
# Add CPPFLAGS to support use of QHull for E-field interpolation

name := G4cmp

ifndef G4INSTALL
$(error Please configure your GEANT4 environment before building G4CMP)
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -Dqh_QHpointer -I$(G4WORKDIR)/include -I$(G4INCLUDE)

# Library goes into user's G4WORKDIR, not into system install
G4LIBDIR := $(G4WORKDIR)/lib/$(G4SYSTEM)

include $(G4INSTALL)/config/common.gmk
