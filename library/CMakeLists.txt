#----------------------------------------------------------------------------
# Need to include G4CMP directories before G4.
#
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

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

#----------------------------------------------------------------------------
# Qhull Compilation stuff
#
set(
    libqhull_HEADERS
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhull/libqhull.h
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhull/geom.h
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhull/io.h
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhull/mem.h
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhull/merge.h
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhull/poly.h
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhull/qhull_a.h
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhull/qset.h
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhull/random.h
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhull/stat.h
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhull/user.h
)
set(
    libqhull_SOURCES
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhull/global.c
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhull/stat.c
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhull/geom2.c
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhull/poly2.c
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhull/merge.c
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhull/libqhull.c
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhull/geom.c
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhull/poly.c
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhull/qset.c
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhull/mem.c
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhull/random.c
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhull/usermem.c
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhull/userprintf.c
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhull/io.c
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhull/user.c
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhull/rboxlib.c
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhull/userprintf_rbox.c
        ${libqhull_HEADERS}
)
set(
    libqhullcpp_HEADERS
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/Coordinates.h
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/functionObjects.h
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/PointCoordinates.h
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/Qhull.h
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullError.h
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullFacet.h
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullFacetList.h
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullFacetSet.h
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullHyperplane.h
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullIterator.h
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullLinkedList.h
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullPoint.h
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullPoints.h
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullPointSet.h
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullQh.h
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullRidge.h
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullSet.h
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullSets.h
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullStat.h
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullVertex.h
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullVertexSet.h
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/RboxPoints.h
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/RoadError.h
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/RoadLogEvent.h
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/UsingLibQhull.h
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/qhulltest/RoadTest.h
)

set(
    libqhullcpp_SOURCES
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/Coordinates.cpp
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/PointCoordinates.cpp
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/Qhull.cpp
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullFacet.cpp
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullFacetList.cpp
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullFacetSet.cpp
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullHyperplane.cpp
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullPoint.cpp
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullPointSet.cpp
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullPoints.cpp
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullQh.cpp
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullRidge.cpp
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullSet.cpp
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullStat.cpp
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullVertex.cpp
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullVertexSet.cpp
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/RboxPoints.cpp
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/RoadError.cpp
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/RoadLogEvent.cpp
        ${PROJECT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/UsingLibQhull.cpp
        ${libqhullcpp_HEADERS}
)
add_library(qhullstatic_p SHARED ${libqhull_SOURCES})
set_target_properties(qhullstatic_p PROPERTIES
    COMPILE_DEFINITIONS "qh_QHpointer"
    VERSION 6.3.1.1494)

add_library(qhullcpp SHARED ${libqhullcpp_SOURCES})
set_target_properties(qhullcpp PROPERTIES
    COMPILE_DEFINITIONS "qh_QHpointer"
    VERSION 6.3.1.1494)
    
include_directories(${PROJECT_SOURCE_DIR}/qhull-2012.1/src)

set(library_SOURCES 
   ${PROJECT_SOURCE_DIR}/library/src/G4CMPDriftBoundaryProcess.cc
   ${PROJECT_SOURCE_DIR}/library/src/G4CMPDriftElectron.cc
   ${PROJECT_SOURCE_DIR}/library/src/G4CMPDriftHole.cc
   ${PROJECT_SOURCE_DIR}/library/src/G4CMPElectrodeHit.cc
   ${PROJECT_SOURCE_DIR}/library/src/G4CMPElectrodeSensitivity.cc
   ${PROJECT_SOURCE_DIR}/library/src/G4CMPeLukeScattering.cc
   ${PROJECT_SOURCE_DIR}/library/src/G4CMPEqEMField.cc
   ${PROJECT_SOURCE_DIR}/library/src/G4CMPFieldManager.cc
   ${PROJECT_SOURCE_DIR}/library/src/G4CMPhLukeScattering.cc
   ${PROJECT_SOURCE_DIR}/library/src/G4CMPInterValleyScattering.cc
   ${PROJECT_SOURCE_DIR}/library/src/G4CMPLocalElectroMagField.cc
   ${PROJECT_SOURCE_DIR}/library/src/G4CMPPhysics.cc
   ${PROJECT_SOURCE_DIR}/library/src/G4CMPPhysicsList.cc
   ${PROJECT_SOURCE_DIR}/library/src/G4CMPProcessUtils.cc
   ${PROJECT_SOURCE_DIR}/library/src/G4CMPStackingAction.cc
   ${PROJECT_SOURCE_DIR}/library/src/G4CMPTimeStepper.cc
   ${PROJECT_SOURCE_DIR}/library/src/G4CMPTriLinearInterp.cc
   ${PROJECT_SOURCE_DIR}/library/src/G4CMPValleyTrackMap.cc
   ${PROJECT_SOURCE_DIR}/library/src/G4CMPVDriftProcess.cc
   ${PROJECT_SOURCE_DIR}/library/src/G4LatticeLogical.cc
   ${PROJECT_SOURCE_DIR}/library/src/G4LatticeManager.cc
   ${PROJECT_SOURCE_DIR}/library/src/G4LatticePhysical.cc
   ${PROJECT_SOURCE_DIR}/library/src/G4LatticeReader.cc
   ${PROJECT_SOURCE_DIR}/library/src/G4PhononDownconversion.cc
   ${PROJECT_SOURCE_DIR}/library/src/G4PhononLong.cc
   ${PROJECT_SOURCE_DIR}/library/src/G4PhononPolarization.cc
   ${PROJECT_SOURCE_DIR}/library/src/G4PhononReflection.cc
   ${PROJECT_SOURCE_DIR}/library/src/G4PhononScattering.cc
   ${PROJECT_SOURCE_DIR}/library/src/G4PhononTrackMap.cc
   ${PROJECT_SOURCE_DIR}/library/src/G4PhononTransFast.cc
   ${PROJECT_SOURCE_DIR}/library/src/G4PhononTransSlow.cc
   ${PROJECT_SOURCE_DIR}/library/src/G4VPhononProcess.cc
   )
    
add_library(G4cmp SHARED ${library_SOURCES})
target_link_libraries(G4cmp ${Geant4_LIBRARIES} qhullcpp qhullstatic_p)
set_target_properties(G4cmp PROPERTIES
                        COMPILE_DEFINITIONS "qh_QHpointer")
set(CMAKE_CXX_FLAGS "-O2 -p -g")
install(TARGETS G4cmp DESTINATION lib)