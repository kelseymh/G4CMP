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
#
# NOTE: Tests do not need RPATH modification because they don't get installed.
add_executable(electron_Epv electron_Epv.cc)
target_link_libraries(electron_Epv G4cmp)

add_executable(latticeVecs latticeVecs.cc)
target_link_libraries(latticeVecs G4cmp)

add_executable(luke_dist luke_dist.cc)
target_link_libraries(luke_dist G4cmp)

add_executable(testBlockData testBlockData.cc)
target_link_libraries(testBlockData G4cmp)

add_executable(testCrystalGroup testCrystalGroup.cc)
target_link_libraries(testCrystalGroup G4cmp)

add_executable(g4cmpEFieldTest g4cmpEFieldTest.cc)
target_link_libraries(g4cmpEFieldTest G4cmp)

add_executable(phononKinematics phononKinematics.cc)
target_link_libraries(phononKinematics G4cmp)

add_executable(testChargeCloud testChargeCloud.cc)
target_link_libraries(testChargeCloud G4cmp)

add_executable(testPartition testPartition.cc)
target_link_libraries(testPartition G4cmp)

add_executable(testHVtransform testHVtransform.cc)
target_link_libraries(testHVtransform G4cmp)
