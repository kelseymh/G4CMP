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
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/library/include)
include(${Geant4_USE_FILE})

set(phonon_SOURCES 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/PhononActionInitialization.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/PhononDetectorConstruction.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/PhononPrimaryGeneratorAction.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/PhononSensitivity.cc
   )
    
add_executable(g4cmpPhonon g4cmpPhonon.cc ${phonon_SOURCES})
target_link_libraries(g4cmpPhonon G4cmp)
set_target_properties(g4cmpPhonon PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib"
    INSTALL_RPATH_USE_LINK_RPATH TRUE
    BUILD_WITH_INSTALL_RPATH FALSE)
install(TARGETS g4cmpPhonon DESTINATION bin)
