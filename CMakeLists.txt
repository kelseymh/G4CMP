#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(G4CMP)

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
    add_subdirectory(channeling)
endif()
