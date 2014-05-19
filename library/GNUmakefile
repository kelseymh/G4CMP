# G4CMP/library/GNUmakefile	Build common-code library libG4cmp.(a|so|dylib)
# $Id$

# Follow GEANT4 library building convention
# Add CPPFLAGS to support use of QHull for E-field interpolation
# Add G4CMP_DEBUG support
# Add "link-time" dependendence on QHull

name := G4cmp

ifndef G4INSTALL
$(error Please configure your GEANT4 environment before building G4CMP)
endif

include $(G4INSTALL)/config/architecture.gmk

G4CMP_FLAGS += -I$(G4INCLUDE) -I$(G4WORKDIR)/include -Dqh_QHpointer 
ifdef G4CMP_DEBUG
  G4CMP_FLAGS += -DG4CMP_DEBUG
endif

CPPFLAGS := $(G4CMP_FLAGS) $(CPPFLAGS)		# Prepend for right load order

# Library goes into user's G4WORKDIR, not into system install
G4LIBDIR := $(G4WORKDIR)/lib/$(G4SYSTEM)

include $(G4INSTALL)/config/common.gmk

# Include linking against Qhull when forming G4CMP shared library
# FIXME:  Needed on MacOSX 10.5.8 (GCC 4.0.1), not other platforms
G4CMP_LIBDEP := -L$(G4WORKDIR)/lib/$(G4SYSTEM) -lqhullcpp -lqhullstatic_p
G4CMP_LIBDEP += -Wl,-rpath,$(G4WORKDIR)/lib/$(G4SYSTEM)
INTYLIBS += $(G4CMP_LIBDEP)
