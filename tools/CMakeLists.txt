#----------------------------------------------------------------------------
# Need to include G4CMP directories before G4.
#
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
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
# Executables include supporting source files from include/ and src/
#

set(kvmap_SOURCES 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPInterpolator.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPPhononKVgTable.cc
    )

add_executable(g4cmpKVtables g4cmpKVtables.cc $(kvmap_SOURCES))
target_link_libraries(g4cmpKVtables G4cmp)
set_target_properties(g4cmpKVtables PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib"
    INSTALL_RPATH_USE_LINK_RPATH TRUE
    BUILD_WITH_INSTALL_RPATH FALSE)

