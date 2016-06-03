#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(G4CMP)

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
# Options to build examples
#

option(BUILD_ALL_EXAMPLES "Build phonon, charge, and charge channeling
examples. Default: OFF" OFF)
option(BUILD_PHONON_EXAMPLE "Build phonon example. Default: OFF" OFF)
option(BUILD_CHARGE_EXAMPLE "Build charge example. Default: OFF" OFF)
option(BUILD_CHANNELING_EXAMPLE "Build charge channeling example. Default: OFF" OFF)

#-----------------------------------------------------------------------------
# Add subdirectories based on options selected
#

add_subdirectory(library)

if(BUILD_ALL_EXAMPLES OR BUILD_PHONON_EXAMPLE)
    add_subdirectory(phonon)
endif()

if(BUILD_ALL_EXAMPLES OR BUILD_CHARGE_EXAMPLE)
    add_subdirectory(charge)
endif()

if(BUILD_ALL_EXAMPLES OR BUILD_CHANNELING_EXAMPLE)
    #FIXME: Channeling demo is broken right now
    #add_subdirectory(channeling)
endif()

#-----------------------------------------------------------------------------
# Install CMake config files
#

install(FILES ${PROJECT_SOURCE_DIR}/cmake/G4CMPConfig.cmake
        DESTINATION cmake)
install(FILES ${PROJECT_SOURCE_DIR}/cmake/UseG4CMP.cmake
        DESTINATION cmake)
