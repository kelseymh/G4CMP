#----------------------------------------------------------------------------
# Find Geant4 package
# NOTE: WITH_GEANT4_UIVIS and USE_GEANT4_STATIC_LIBS are defined here
#
if(NOT Geant4_FOUND)
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

set(sensor_SOURCES 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ChargeFETDigitizerModule.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ChargeFETDigitizerMessenger.cc
    )

set(fet_CONFIGS
    ${CMAKE_CURRENT_SOURCE_DIR}/config/ConstantsFET
    ${CMAKE_CURRENT_SOURCE_DIR}/config/FETTemplates
    ${CMAKE_CURRENT_SOURCE_DIR}/config/EpotRamoChan1
    ${CMAKE_CURRENT_SOURCE_DIR}/config/EpotRamoChan2
    ${CMAKE_CURRENT_SOURCE_DIR}/config/EpotRamoChan3
    ${CMAKE_CURRENT_SOURCE_DIR}/config/EpotRamoChan4
    )

if(USE_GEANT4_STATIC_LIBS)
    add_library(G4cmpSensor STATIC ${sensor_SOURCES})
else()
    add_library(sensorLib SHARED ${sensor_SOURCES})
endif()
set_target_properties(sensorLib PROPERTIES OUTPUT_NAME G4cmpSensor)
    
target_link_libraries(sensorLib G4cmp)

add_executable(g4cmpFETSim g4cmpFETSim.cc)
target_link_libraries(g4cmpFETSim sensorLib)

install(TARGETS sensorLib DESTINATION lib)
install(TARGETS g4cmpFETSim DESTINATION bin)
install(FILES ${fet_CONFIGS} DESTINATION config/G4CMP/FETSim)
