#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(G4CMP)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Options to build examples
#

option(BUILD_ALL_EXAMPLES "Build phonon, charge, and charge channeling
examples. Default: OFF" OFF)
option(BUILD_PHONON_EXAMPLE "Build phonon example. Default: OFF" OFF)
option(BUILD_CHARGE_EXAMPLE "Build charge example. Default: OFF" OFF)
option(BUILD_CHANNELING_EXAMPLE "Build charge channeling example. Default: OFF" OFF)
#option(BUILD_LIBRARY "Build G4CMP LIBRARY. Default: ON" ON)


#----------------------------------------------------------------------------
# Grab Qhull files needed to link against it.
# 

include_directories(${PROJECT_SOURCE_DIR}/qhull-2012.1/src)


#--------------------------------------------------------------------------------------
# Add executables and link to Geant4 and Qhull
#

add_subdirectory(library)

if(BUILD_ALL_EXAMPLES OR BUILD_PHONON_EXAMPLE)
    add_subdirectory(phonon)
endif()

if(BUILD_ALL_EXAMPLES OR BUILD_CHARGE_EXAMPLE)
    add_subdirectory(charge)
endif()

if(BUILD_ALL_EXAMPLES OR BUILD_CHANNELING_EXAMPLE)
    add_subdirectory(channeling)
endif()
