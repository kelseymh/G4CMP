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

#----------------------------------------------------------------------------
# RPATH stuff
#
# Automatically add paths of linked libraries to RPATH *unless* they are in the 
# build directory
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# Set the RPATH of the installed objects *unless* they are installing to a 
# system directory, in which case they shouldn't need an explicitly set RPATH.
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES 
     "${CMAKE_INSTALL_PREFIX}/lib" is_SYSTEM_DIR)
if("${is_SYSTEM_DIR}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif()

#----------------------------------------------------------------------------
# Build stuff
#
set(phonon_SOURCES 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/PhononActionInitialization.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/PhononDetectorConstruction.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/PhononPrimaryGeneratorAction.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/PhononSensitivity.cc
    )
    
add_library(G4cmpPhonon SHARED ${phonon_SOURCES})
target_link_libraries(G4cmpPhonon G4cmp)

add_executable(g4cmpPhonon g4cmpPhonon.cc)
target_link_libraries(g4cmpPhonon G4cmpPhonon)

install(TARGETS G4cmpPhonon DESTINATION lib)
install(TARGETS g4cmpPhonon DESTINATION bin)
