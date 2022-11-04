# G4CMP/library/GNUmakefile	Build common-code library libG4cmp.(a|so|dylib)
# $Id: a3356a2caa936d64b7c9ed2c6112d8e18aca3a6b $

# Follow GEANT4 library building convention
# Add CPPFLAGS to support use of Qhull for E-field interpolation
# Add G4CMP_DEBUG support
# Add "link-time" dependendence on Qhull
# Add G4CMP_SET_ELECTRON_MASS support to build with dynamic mass code
# Add qhull/src directory to include directories list
# Use G4CMPINSTALL instead of ".." to find includes
# Move QHull building under library
# Drop G4CMP_SET_ELECTRON_MASS code blocks; not physical
# Ensure that G4LIBDIR is set correctly, overriding common.gmk.
# Add -DG4DIGI_ALLOC_EXPORT flag to support new hits collection
# Add G4CMP_USE_SANITIZER, G4CMP_SANITIZER_TYPE for thread-safety checking
# Add G4LIB_USE_CLHEP to distinguish G4's DoubConv.h from CLHEP's DoubConv.hh
# Use G4DEBUG to select optimization level; include debugging symbols always

name := G4cmp

ifndef G4INSTALL
$(error Please configure your GEANT4 environment before building G4CMP)
endif

include $(G4INSTALL)/config/architecture.gmk

G4CMP_FLAGS += -I./include -I$(G4CMPINSTALL)/library/qhull-2012.1/src \
	-I$(G4INCLUDE) -Dqh_QHpointer -DG4DIGI_ALLOC_EXPORT

ifdef G4CMP_DEBUG
  G4CMP_FLAGS += -DG4CMP_DEBUG
endif
ifdef G4CMPTLI_DEBUG
  G4CMP_FLAGS += -DG4CMPTLI_DEBUG
endif
ifdef G4CMP_USE_SANITIZER
  G4CMP_SANITIZER_TYPE := thread		# User can override w/envvar
  G4CMP_FLAGS += -fno-omit-frame-pointer -fsanitize=$(G4CMP_SANITIZER_TYPE)
endif
ifdef G4LIB_USE_CLHEP				# Using internal CLHEP sources
  G4CMP_FLAGS += -DG4LIB_USE_CLHEP
endif

CPPFLAGS := $(G4CMP_FLAGS) $(CPPFLAGS)		# Prepend for right load order
ifdef G4DEBUG
  CXXFLAGS += -Og -g -pipe	# ??What does the "-pipe" flag actually do??
else
  CXXFLAGS += -O2 -g -pipe
endif

# Library goes into user's G4WORKDIR, not into system install
override G4LIBDIR := $(G4WORKDIR)/lib/$(G4SYSTEM)

.PHONY : lib clean

lib : qhull		# Set dependence on QHull before Geant4 configuration
clean :: qhull.cleanall

include $(G4INSTALL)/config/common.gmk

# Include linking against Qhull when forming G4CMP shared library
# FIXME:  Needed on MacOSX 10.5.8 (GCC 4.0.1), not other platforms
G4CMP_LIBDEP := -L$(G4LIBDIR) -lqhullcpp -lqhullstatic_p
G4CMP_LIBDEP += -Wl,-rpath,$(G4LIBDIR)
INTYLIBS += $(G4CMP_LIBDEP)

# Manually configure building the Qhull libraries in Geant4 style
ISMAC := $(findstring Darwin,$(G4SYSTEM))
ISWIN := $(findstring Win,$(G4SYSTEM))
DYLIB_OPTS := -dynamiclib -undefined suppress -flat_namespace

.PHONY : qhull qhull.all qhull.install qhull.cleanall

qhull :
	-$(MAKE) -C qhull-2012.1 DESTDIR=$(G4WORKDIR) \
	  LIBDIR=$(G4LIBDIR) \
	  BINDIR=$(G4WORKDIR)/bin/$(G4SYSTEM) \
	  CC_OPTS3="$(if $(ISMAC),$(DYLIB_OPTS),)" \
	  SO=$(if $(ISMAC),dylib,$(if $(ISWIN),dll,so)).6.3.1 \
	  all install

qhull.% :
	-$(MAKE) -C qhull-2012.1 $(subst .,,$(suffix $@))
