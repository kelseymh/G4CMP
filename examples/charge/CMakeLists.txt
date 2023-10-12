cmake_minimum_required(VERSION 3.0)
project(g4cmpCharge)

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

set(charge_SOURCES 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ChargeActionInitialization.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ChargeConfigManager.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ChargeConfigMessenger.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ChargeDetectorConstruction.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ChargeElectrodePattern.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ChargeElectrodeSensitivity.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ChargePrimaryGeneratorAction.cc
    )

set(charge_EPOT_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/EPotFiles/EPot_iZIP4_4V
    ${CMAKE_CURRENT_SOURCE_DIR}/EPotFiles/EPot_iZIP4_4V_small
    ${CMAKE_CURRENT_SOURCE_DIR}/EPotFiles/EPot_iZIP4_70V_small
    ${CMAKE_CURRENT_SOURCE_DIR}/EPotFiles/EPot_iZIP4_axissymmetric_shortQOS1
    )

set(charge_MACROS
    ${CMAKE_CURRENT_SOURCE_DIR}/run.mac
    ${CMAKE_CURRENT_SOURCE_DIR}/vis.mac
    ${CMAKE_CURRENT_SOURCE_DIR}/fast.mac
    ${CMAKE_CURRENT_SOURCE_DIR}/loop.mac
    ${CMAKE_CURRENT_SOURCE_DIR}/open.mac
    ${CMAKE_CURRENT_SOURCE_DIR}/play.mac
    ${CMAKE_CURRENT_SOURCE_DIR}/iZip4.mac
    ${CMAKE_CURRENT_SOURCE_DIR}/single.mac
    ${CMAKE_CURRENT_SOURCE_DIR}/openNew.mac
    ${CMAKE_CURRENT_SOURCE_DIR}/useRndm.mac
    ${CMAKE_CURRENT_SOURCE_DIR}/saveRndm.mac
    ${CMAKE_CURRENT_SOURCE_DIR}/testRefl.mac
    ${CMAKE_CURRENT_SOURCE_DIR}/shortQOS1.mac
    ${CMAKE_CURRENT_SOURCE_DIR}/testIonize.mac
    ${CMAKE_CURRENT_SOURCE_DIR}/drift_curve.mac
    ${CMAKE_CURRENT_SOURCE_DIR}/movie.loop
    )

if(USE_GEANT4_STATIC_LIBS)
    add_library(chargeLib STATIC ${charge_SOURCES})
else()
    add_library(chargeLib SHARED ${charge_SOURCES})
endif()
set_target_properties(chargeLib PROPERTIES OUTPUT_NAME g4cmpCharge)
    
target_link_libraries(chargeLib  ${G4CMP_LIBRARIES} ${Geant4_LIBRARIES})

add_executable(g4cmpCharge g4cmpCharge.cc)
target_link_libraries(g4cmpCharge chargeLib)

add_executable(g4cmpIonize g4cmpIonize.cc)
target_link_libraries(g4cmpIonize chargeLib)

install(TARGETS chargeLib DESTINATION lib)
install(TARGETS g4cmpCharge DESTINATION bin)
install(TARGETS g4cmpIonize DESTINATION bin)
install(FILES ${charge_EPOT_FILES} DESTINATION EPotFiles COMPONENT config)
install(FILES ${charge_MACROS} DESTINATION macros/g4cmpCharge)
