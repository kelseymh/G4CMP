#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(G4CMP)

# User can provide version string in case build is not from Git source
set(G4CMP_VERSION "" CACHE STRING "G4CMP Version Tag to be used")

#----------------------------------------------------------------------------
# Find Geant4 package
# NOTE: WITH_GEANT4_UIVIS and USE_GEANT4_STATIC_LIBS are defined here
#
include(FindGeant4.cmake)

#----------------------------------------------------------------------------
# Set up compiler flags for various build types
#
# Always need to build against C++11
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS_DEV "${CMAKE_CXX_FLAGS} -g -p -O2 -Wall -Wextra -Wshadow -Wpedantic")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG} -p -DG4CMP_DEBUG=1")
set(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_MINSIZEREL}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE} -pipe")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -pipe -p")

# Set default build type as "Release"
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

#----------------------------------------------------------------------------
# Utility function to define simple executables in subdirectories
#
function(make_binaries)
    foreach(bin IN LISTS ARGV)
        add_executable(${bin} ${bin}.cc)
        target_link_libraries(${bin} G4cmp)
	install(TARGETS ${bin} DESTINATION bin)
    endforeach()
endfunction()

#----------------------------------------------------------------------------
# Options to build examples and other directories
#
option(BUILD_G4CMP_TOOLS "Build utility and support programs.  Default: ON" ON)
option(BUILD_G4CMP_TESTS "Build unit tests for classes.  Default: OFF" OFF)
option(BUILD_ALL_EXAMPLES "Build phonon, charge and sensor examples. Default: OFF" OFF)
option(BUILD_PHONON_EXAMPLE "Build phonon example. Default: OFF" OFF)
option(BUILD_CHARGE_EXAMPLE "Build charge example. Default: OFF" OFF)
option(BUILD_SENSOR_EXAMPLE "Build sensor example. Default: OFF" OFF)

#-----------------------------------------------------------------------------
# Add subdirectories based on options selected
#
# Always build the library
add_subdirectory(library)

if (BUILD_G4CMP_TOOLS)
    add_subdirectory(tools)
endif()

if (BUILD_G4CMP_TESTS)
    add_subdirectory(tests)
endif()

#-----------------------------------------------------------------------------
# Examples: built using extra 'make' targets; flags include with 'all'
#

if(0)	# TEMPORARILY DISABLE EXAMPLES
add_custom_target(examples $<BUILD_ALL_EXAMPLES:ALL>
    DEPENDS phonon charge sensor)

add_custom_target(phonon $<$<OR:BUILD_ALL_EXAMPLES,BUILD_PHONON_EXAMPLE>:ALL>)
add_subdirectory(examples/phonon EXCLUDE_FROM_ALL)
add_dependencies(phonon install examples/phonon)

add_custom_target(charge $<$<OR:BUILD_ALL_EXAMPLES,BUILD_CHARGE_EXAMPLE>:ALL>)
add_subdirectory(examples/charge EXCLUDE_FROM_ALL)
add_dependencies(charge install examples/charge)

add_custom_target(sensor $<$<OR:BUILD_ALL_EXAMPLES,BUILD_SENSOR_EXAMPLE>:ALL>)
add_subdirectory(examples/sensor EXCLUDE_FROM_ALL)
add_dependencies(sensor install examples/sensor)
endif()

#-----------------------------------------------------------------------------
# Create a version file as part of the "make all" procedure

set(version_FILE .g4cmp-version)

find_program(SH sh HINTS /bin)
add_custom_target(version ALL
    COMMENT Building version identification file
    VERBATIM			# Allows for parentheses and redirection below
    COMMAND ${CMAKE_COMMAND} -E env ${SH}
        -c "(cd ${PROJECT_SOURCE_DIR};\
        [ -d .git ] && git describe || echo ${G4CMP_VERSION}) > ${version_FILE}"
    )

install(FILES ${PROJECT_BINARY_DIR}/${version_FILE} DESTINATION share/G4CMP COMPONENT config)

#-----------------------------------------------------------------------------
# Install CMake config files
#
install(FILES ${PROJECT_SOURCE_DIR}/cmake/G4CMPConfig.cmake
        DESTINATION cmake)
install(FILES ${PROJECT_SOURCE_DIR}/cmake/UseG4CMP.cmake
        DESTINATION cmake)
