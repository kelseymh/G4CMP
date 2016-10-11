#----------------------------------------------------------------------------
# Find Geant4 package
# NOTE: WITH_GEANT4_UIVIS and USE_GEANT4_STATIC_LIBS are defined here
#
if(NOT ${Geant4_FOUND})
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
add_executable(electron_Epv electron_Epv.cc)
target_link_libraries(electron_Epv G4cmp)
set_target_properties(electron_Epv PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib"
    INSTALL_RPATH_USE_LINK_RPATH TRUE
    BUILD_WITH_INSTALL_RPATH FALSE)

add_executable(latticeVecs latticeVecs.cc)
target_link_libraries(latticeVecs G4cmp)
set_target_properties(latticeVecs PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib"
    INSTALL_RPATH_USE_LINK_RPATH TRUE
    BUILD_WITH_INSTALL_RPATH FALSE)

add_executable(luke_dist luke_dist.cc)
target_link_libraries(luke_dist G4cmp)
set_target_properties(luke_dist PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib"
    INSTALL_RPATH_USE_LINK_RPATH TRUE
    BUILD_WITH_INSTALL_RPATH FALSE)

add_executable(testBlockData testBlockData.cc)
target_link_libraries(testBlockData G4cmp)
set_target_properties(testBlockData PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib"
    INSTALL_RPATH_USE_LINK_RPATH TRUE
    BUILD_WITH_INSTALL_RPATH FALSE)

add_executable(testCrystalGroup testCrystalGroup.cc)
target_link_libraries(testCrystalGroup G4cmp)
set_target_properties(testCrystalGroup PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib"
    INSTALL_RPATH_USE_LINK_RPATH TRUE
    BUILD_WITH_INSTALL_RPATH FALSE)

add_executable(g4cmpEFieldTest g4cmpEFieldTest.cc)
target_link_libraries(g4cmpEFieldTest G4cmp)
set_target_properties(g4cmpEFieldTest PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib"
    INSTALL_RPATH_USE_LINK_RPATH TRUE
    BUILD_WITH_INSTALL_RPATH FALSE)
