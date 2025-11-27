cmake_minimum_required(VERSION 3.0)
project(g4cmpSensor)

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

# Use same C++ version as was used for Geant4 (v10 or v11)
if (NOT DEFINED CMAKE_CXX_STANDARD)	# Only if user did not specify
    if (DEFINED Geant4_CXXSTD)			# For G4 v10, value is "c++NN"
        string(REPLACE "c++" "" CMAKE_CXX_STANDARD ${Geant4_CXXSTD})
    elseif (DEFINED Geant4_CXX_STANDARD)	# For G4 v11, value is "NN"
        set(CMAKE_CXX_STANDARD ${Geant4_CXX_STANDARD})
    endif()
endif()

set(CMAKE_CXX_STANDARD_REQUIRED True)

# Store actual value of c++NN derived from above or set by user
set(G4CMP_CXX_STANDARD ${CMAKE_CXX_STANDARD} CACHE STRING
   "C++ version used for building G4CMP")

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
    
target_link_libraries(sensorLib ${G4CMP_LIBRARIES} ${Geant4_LIBRARIES})

add_executable(g4cmpFETSim g4cmpFETSim.cc)
target_link_libraries(g4cmpFETSim sensorLib)

install(TARGETS sensorLib DESTINATION lib)
install(TARGETS g4cmpFETSim DESTINATION bin)
install(FILES ${fet_CONFIGS} DESTINATION config/G4CMP/FETSim)
