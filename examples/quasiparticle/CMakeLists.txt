cmake_minimum_required(VERSION 3.0)
project(g4cmpQuasiparticle)

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
set(quasiparticle_SOURCES 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/QuasiparticleActionInitialization.cc 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/QuasiparticleConfigManager.cc 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/QuasiparticleConfigMessenger.cc 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/QuasiparticleDetectorConstruction.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/QuasiparticlePrimaryGeneratorAction.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/QuasiparticleSensitivity.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/QuasiparticleSteppingAction.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/QuasiparticleQubitHousing.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/QuasiparticlePad.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/QuasiparticleTransmissionLine.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/QuasiparticleResonatorAssembly.cc
    )
    
if(USE_GEANT4_STATIC_LIBS)
    add_library(quasiparticleLib STATIC ${quasiparticle_SOURCES})
else()
    add_library(quasiparticleLib SHARED ${quasiparticle_SOURCES})
endif()
set_target_properties(quasiparticleLib PROPERTIES OUTPUT_NAME g4cmpQuasiparticle)

target_link_libraries(quasiparticleLib ${G4CMP_LIBRARIES} ${Geant4_LIBRARIES})

add_executable(g4cmpQuasiparticle g4cmpQuasiparticle.cc)
target_link_libraries(g4cmpQuasiparticle quasiparticleLib)

install(TARGETS quasiparticleLib DESTINATION lib)
install(TARGETS g4cmpQuasiparticle DESTINATION bin)
