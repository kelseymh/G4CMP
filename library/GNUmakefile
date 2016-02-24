# G4CMP/library/GNUmakefile	Build common-code library libG4cmp.(a|so|dylib)
# $Id: a3356a2caa936d64b7c9ed2c6112d8e18aca3a6b $

# Follow GEANT4 library building convention
# Add CPPFLAGS to support use of Qhull for E-field interpolation
# Add G4CMP_DEBUG support
# Add "link-time" dependendence on Qhull
# Add G4CMP_SET_ELECTRON_MASS support to build with dynamic mass code

name := G4cmp

ifndef G4INSTALL
$(error Please configure your GEANT4 environment before building G4CMP)
endif

include $(G4INSTALL)/config/architecture.gmk

G4CMP_FLAGS += -I./include -I$(G4WORKDIR)/include -I$(G4INCLUDE) -Dqh_QHpointer 
ifdef G4CMP_DEBUG
  G4CMP_FLAGS += -DG4CMP_DEBUG
endif
ifdef G4CMP_SET_ELECTRON_MASS
  G4CMP_FLAGS += -DG4CMP_SET_ELECTRON_MASS
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
