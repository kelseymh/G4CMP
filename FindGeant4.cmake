#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#

# NOTE: Geant4's CMake file treats ALL components as REQUIRED if we set
# the package as REQUIRED. In other words, OPTIONAL_COMPONENTS are still
# seen as REQUIRED and CMake will fail if we try to find optional components.
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
option(USE_GEANT4_STATIC_LIBS "Build example with Geant4 static libraries" OFF)

if(WITH_GEANT4_UIVIS)
    list(APPEND COMPONENTS ui_all vis_all)
endif()

if(USE_GEANT4_STATIC_LIBS)
    list(APPEND COMPONENTS static)
endif()

# G4CMP Requires G4 version 10.4 or greater.
find_package(Geant4 10.4 REQUIRED ${COMPONENTS})

set(Geant4_FOUND TRUE CACHE BOOL "FindGeant4.cmake was successful")
