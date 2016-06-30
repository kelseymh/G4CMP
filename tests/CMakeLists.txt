#----------------------------------------------------------------------------
# Need to include G4CMP directories before G4.
#
include_directories(${PROJECT_SOURCE_DIR}/library/include)

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

add_executable(testG4CMPBlockData testG4CMPBlockData.cc)
target_link_libraries(testBlockData G4cmp)
set_target_properties(testBlockData PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib"
    INSTALL_RPATH_USE_LINK_RPATH TRUE
    BUILD_WITH_INSTALL_RPATH FALSE)
