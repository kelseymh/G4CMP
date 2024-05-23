#----------------------------------------------------------------------------
# Find Geant4 package
# NOTE: WITH_GEANT4_UIVIS and USE_GEANT4_STATIC_LIBS are defined here
#
if(NOT Geant4_FOUND)
    include(${PROJECT_SOURCE_DIR}/FindGeant4.cmake)
endif()

#----------------------------------------------------------------------------
# Setup include directories and compile definitions
# NOTE: Need to include G4CMP directories before G4.
#
include_directories(${PROJECT_SOURCE_DIR}/library/include)
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Executables are single-file builds, with no associated local library
# NOTE: Add names of binaries to list
#
make_binaries("electron_Epv" "latticeVecs" "luke_dist" "testBlockData"
              "testCrystalGroup" "g4cmpEFieldTest"
              "testChargeCloud" "testPartition" "testHVtransform"
              "testFanoFactor" "testTemperature" )

