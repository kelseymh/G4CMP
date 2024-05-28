cmake_minimum_required(VERSION 3.0)
project(RISQTutorial)

#----------------------------------------------------------------------------
# Find G4CMP package -- NOTE: Assumes user has run g4cmp_env.*sh script!
#
set(CMAKE_PREFIX_PATH $ENV{G4CMPINSTALL}/../..)
find_package(G4CMP REQUIRED)
include(${G4CMP_USE_FILE})

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
set(RISQTutorial_SOURCES 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/RISQTutorialActionInitialization.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/RISQTutorialSteppingAction.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/RISQTutorialConfigManager.cc 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/RISQTutorialConfigMessenger.cc 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/RISQTutorialDetectorConstruction.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/RISQTutorialPrimaryGeneratorAction.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/RISQTutorialSensitivity.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/RISQTutorialQubitHousing.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/RISQTutorialPad.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/RISQTutorialTransmissionLine.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/RISQTutorialResonatorAssembly.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/RISQTutorialStraightFluxLine.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/RISQTutorialCornerFluxLine.cc
    )
    
if(USE_GEANT4_STATIC_LIBS)
    add_library(RISQTutorialLib STATIC ${RISQTutorial_SOURCES})
else()
    add_library(RISQTutorialLib SHARED ${RISQTutorial_SOURCES})
endif()
set_target_properties(RISQTutorialLib PROPERTIES OUTPUT_NAME g4cmpRISQTutorial)

message("G4CMP Libraries: ")
message(${G4CMP_LIBRARIES})

target_link_libraries(RISQTutorialLib ${G4CMP_LIBRARIES} ${Geant4_LIBRARIES})

add_executable(RISQTutorial RISQTutorial.cc)
target_link_libraries(RISQTutorial RISQTutorialLib)

install(TARGETS RISQTutorialLib DESTINATION lib)
install(TARGETS RISQTutorial DESTINATION bin)
