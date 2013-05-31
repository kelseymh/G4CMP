#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Qhull Compilation stuff
#
#
set(
    libqhull_HEADERS
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhull/libqhull.h
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhull/geom.h
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhull/io.h
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhull/mem.h
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhull/merge.h
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhull/poly.h
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhull/qhull_a.h
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhull/qset.h
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhull/random.h
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhull/stat.h
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhull/user.h
)
set(
    libqhull_SOURCES
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhull/global.c
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhull/stat.c
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhull/geom2.c
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhull/poly2.c
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhull/merge.c
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhull/libqhull.c
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhull/geom.c
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhull/poly.c
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhull/qset.c
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhull/mem.c
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhull/random.c
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhull/usermem.c
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhull/userprintf.c
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhull/io.c
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhull/user.c
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhull/rboxlib.c
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhull/userprintf_rbox.c
        ${libqhull_HEADERS}
)
set(
    libqhullcpp_HEADERS
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhullcpp/Coordinates.h
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhullcpp/functionObjects.h
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhullcpp/PointCoordinates.h
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhullcpp/Qhull.h
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhullcpp/QhullError.h
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhullcpp/QhullFacet.h
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhullcpp/QhullFacetList.h
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhullcpp/QhullFacetSet.h
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhullcpp/QhullHyperplane.h
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhullcpp/QhullIterator.h
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhullcpp/QhullLinkedList.h
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhullcpp/QhullPoint.h
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhullcpp/QhullPoints.h
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhullcpp/QhullPointSet.h
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhullcpp/QhullQh.h
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhullcpp/QhullRidge.h
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhullcpp/QhullSet.h
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhullcpp/QhullSets.h
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhullcpp/QhullStat.h
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhullcpp/QhullVertex.h
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhullcpp/QhullVertexSet.h
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhullcpp/RboxPoints.h
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhullcpp/RoadError.h
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhullcpp/RoadLogEvent.h
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhullcpp/UsingLibQhull.h
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/qhulltest/RoadTest.h
)

set(
    libqhullcpp_SOURCES
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhullcpp/Coordinates.cpp
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhullcpp/PointCoordinates.cpp
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhullcpp/Qhull.cpp
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhullcpp/QhullFacet.cpp
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhullcpp/QhullFacetList.cpp
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhullcpp/QhullFacetSet.cpp
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhullcpp/QhullHyperplane.cpp
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhullcpp/QhullPoint.cpp
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhullcpp/QhullPointSet.cpp
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhullcpp/QhullPoints.cpp
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhullcpp/QhullQh.cpp
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhullcpp/QhullRidge.cpp
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhullcpp/QhullSet.cpp
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhullcpp/QhullStat.cpp
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhullcpp/QhullVertex.cpp
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhullcpp/QhullVertexSet.cpp
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhullcpp/RboxPoints.cpp
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhullcpp/RoadError.cpp
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhullcpp/RoadLogEvent.cpp
        ${PROJECT_SOURCE_DIR}/include/qhull-2012.1/src/libqhullcpp/UsingLibQhull.cpp
        ${libqhullcpp_HEADERS}
)

add_library(qhullstatic_p STATIC ${libqhull_SOURCES})
set_target_properties(qhullstatic_p PROPERTIES
    COMPILE_DEFINITIONS "qh_QHpointer"
    VERSION 6.3.1.1494)

add_library(qhullcpp STATIC ${libqhullcpp_SOURCES})
set_target_properties(qhullcpp PROPERTIES
    COMPILE_DEFINITIONS "qh_QHpointer"
    VERSION 6.3.1.1494)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(test test.cc ${sources} ${headers})
target_link_libraries(test ${Geant4_LIBRARIES})
target_link_libraries(test qhullcpp)
target_link_libraries(test qhullstatic_p)
set_target_properties(test PROPERTIES
    COMPILE_DEFINITIONS "qh_QHpointer")
set(CMAKE_CXX_FLAGS "-O2")

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS test DESTINATION bin)
