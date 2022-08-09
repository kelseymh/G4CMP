cmake_minimum_required(VERSION 3.0)
project(g4cmpPhonon)

#----------------------------------------------------------------------------
# Find G4CMP package -- NOTE: Assumes user has run g4cmp_env.*sh script!
#
set(CMAKE_PREFIX_PATH $ENV{G4CMPINSTALL}/../..)
find_package(G4CMP)

#----------------------------------------------------------------------------
# Find Geant4 package
# NOTE: WITH_GEANT4_UIVIS and USE_GEANT4_STATIC_LIBS are defined here
#
include($ENV{G4CMPINSTALL}/FindGeant4.cmake)

#----------------------------------------------------------------------------
# Setup include directories and compile definitions
# NOTE: Need to include G4CMP directories before G4.
#
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${G4CMP_INCLUDE_DIRS})
include(${G4CMP_USE_FILE})
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
    ${CMAKE_CURRENT_SOURCE_DIR}/src/PhononConfigManager.cc 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/PhononConfigMessenger.cc 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/PhononDetectorConstruction.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/PhononPrimaryGeneratorAction.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/PhononSensitivity.cc
    )
    
if(USE_GEANT4_STATIC_LIBS)
    add_library(phononLib STATIC ${phonon_SOURCES})
else()
    add_library(phononLib SHARED ${phonon_SOURCES})
endif()
set_target_properties(phononLib PROPERTIES OUTPUT_NAME g4cmpPhonon)

target_link_libraries(phononLib ${G4CMP_LIBRARIES} ${Geant4_LIBRARIES})

add_executable(g4cmpPhonon g4cmpPhonon.cc)
target_link_libraries(g4cmpPhonon phononLib)

install(TARGETS phononLib DESTINATION lib)
install(TARGETS g4cmpPhonon DESTINATION bin)
