cmake_minimum_required(VERSION 3.12)
project(g4cmpValidation)

#----------------------------------------------------------------------------
# Configure include directories and compiler flags for G4CMP and Geant4
# NOTE:  User must source g4cmp_env.sh first, to get envvars

set(CMAKE_PREFIX_PATH $ENV{G4CMPINSTALL}/../..)
find_package(G4CMP)
include(FindGeant4.cmake)

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
set(validation_SOURCES 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ValidationActionInitialization.cc 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ValidationConfigManager.cc 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ValidationConfigMessenger.cc 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ValidationDetectorConstruction.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ValidationPrimaryGeneratorAction.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ValidationSensitivity.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ValidationSteppingAction.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ValidationQubitHousing.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ValidationPad.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ValidationResonatorAssembly.cc	
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ValidationTransmissionLine.cc
    )
    
if(USE_GEANT4_STATIC_LIBS)
    add_library(validationLib STATIC ${validation_SOURCES})
else()
    add_library(validationLib SHARED ${validation_SOURCES})
endif()
set_target_properties(validationLib PROPERTIES OUTPUT_NAME g4cmpValidation)

target_link_libraries(validationLib ${G4CMP_LIBRARIES} ${Geant4_LIBRARIES})

add_executable(g4cmpValidation g4cmpValidation.cc)
target_link_libraries(g4cmpValidation validationLib)

install(TARGETS validationLib DESTINATION lib)
install(TARGETS g4cmpValidation DESTINATION bin)
