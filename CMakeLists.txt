#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(G4CMP)

# User can provide version string in case build is not from Git source
set(G4CMP_VERSION "" CACHE STRING "G4CMP Version Tag to be used")

#----------------------------------------------------------------------------
# Option to include additional diagnostic output
#
option(G4CMP_DEBUG "Include diagnostic output from processes.  Default: OFF" OFF)
if (G4CMP_DEBUG)
    set(CMAKE_CPP_FLAGS "${CMAKE_CPP_FLAGS} -DG4CMP_DEBUG=1")
endif()

#----------------------------------------------------------------------------
# Find Geant4 package
# NOTE: WITH_GEANT4_UIVIS and USE_GEANT4_STATIC_LIBS are defined here
#
include(FindGeant4.cmake)

#----------------------------------------------------------------------------
# Set up compiler flags for various build types
#
# Always need to build against C++11
set(CMAKE_CXX_STANDARD 11)

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
option(INSTALL_EXAMPLES "Copy examples directories to installation area. Default: OFF" OFF)

#-----------------------------------------------------------------------------
# Add subdirectories based on options selected
#
add_subdirectory(library)

if (BUILD_G4CMP_TOOLS)
    add_subdirectory(tools)
endif()

if (BUILD_G4CMP_TESTS)
    add_subdirectory(tests)
endif()

#-----------------------------------------------------------------------------
# Create a version file as part of the "make all" procedure
#
set(version_FILE .g4cmp-version)

find_program(SH sh HINTS /bin)
add_custom_target(version ALL
    COMMENT Building version identification file
    VERBATIM			# Allows for parentheses and redirection below
    COMMAND ${CMAKE_COMMAND} -E env ${SH}
        -c "(cd ${PROJECT_SOURCE_DIR};\
        [ -d .git ] && git describe || echo ${G4CMP_VERSION}) > ${version_FILE}"
    )

install(FILES ${PROJECT_BINARY_DIR}/${version_FILE} DESTINATION share/G4CMP
        COMPONENT config)

#-----------------------------------------------------------------------------
# Install CMake config files
#
install(FILES ${PROJECT_SOURCE_DIR}/cmake/G4CMPConfig.cmake
        DESTINATION cmake)
install(FILES ${PROJECT_SOURCE_DIR}/cmake/UseG4CMP.cmake
        DESTINATION cmake)

#-----------------------------------------------------------------------------
# Install example applications if requested
#
if(INSTALL_EXAMPLES)
    install(DIRECTORY examples DESTINATION . COMPONENT examples
	    PATTERN ".git" EXCLUDE)
endif()
