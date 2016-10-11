#----------------------------------------------------------------------------
# Find Geant4 package
# NOTE: WITH_GEANT4_UIVIS and USE_GEANT4_STATIC_LIBS are defined here
#
if(NOT ${Geant4_FOUND})
    include(${PROJECT_SOURCE_DIR}/FindGeant4.cmake)
endif()

#----------------------------------------------------------------------------
# Setup include directories and compile definitions
# NOTE: Need to include G4CMP directories before G4.
#
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Qhull Compilation stuff -- sources are under current (library/) directory
#
set(libqhull_HEADERS
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhull/libqhull.h
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhull/geom.h
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhull/io.h
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhull/mem.h
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhull/merge.h
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhull/poly.h
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhull/qhull_a.h
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhull/qset.h
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhull/random.h
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhull/stat.h
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhull/user.h
)

set(libqhull_SOURCES
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhull/global.c
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhull/stat.c
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhull/geom2.c
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhull/poly2.c
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhull/merge.c
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhull/libqhull.c
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhull/geom.c
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhull/poly.c
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhull/qset.c
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhull/mem.c
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhull/random.c
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhull/usermem.c
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhull/userprintf.c
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhull/io.c
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhull/user.c
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhull/rboxlib.c
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhull/userprintf_rbox.c
   ${libqhull_HEADERS}
)

set(libqhullcpp_HEADERS
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/Coordinates.h
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/functionObjects.h
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/PointCoordinates.h
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/Qhull.h
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullError.h
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullFacet.h
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullFacetList.h
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullFacetSet.h
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullHyperplane.h
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullIterator.h
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullLinkedList.h
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullPoint.h
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullPoints.h
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullPointSet.h
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullQh.h
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullRidge.h
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullSet.h
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullSets.h
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullStat.h
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullVertex.h
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullVertexSet.h
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/RboxPoints.h
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/RoadError.h
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/RoadLogEvent.h
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/UsingLibQhull.h
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/qhulltest/RoadTest.h
)

set(libqhullcpp_SOURCES
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/Coordinates.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/PointCoordinates.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/Qhull.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullFacet.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullFacetList.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullFacetSet.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullHyperplane.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullPoint.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullPointSet.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullPoints.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullQh.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullRidge.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullSet.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullStat.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullVertex.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/QhullVertexSet.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/RboxPoints.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/RoadError.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/RoadLogEvent.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp/UsingLibQhull.cpp
   ${libqhullcpp_HEADERS}
)

add_library(qhull_p SHARED ${libqhull_SOURCES})
set_target_properties(qhull_p PROPERTIES
    COMPILE_DEFINITIONS "qh_QHpointer"
    VERSION 6.3.1.1494)

add_library(qhullcpp SHARED ${libqhullcpp_SOURCES})
set_target_properties(qhullcpp PROPERTIES
    COMPILE_DEFINITIONS "qh_QHpointer"
    VERSION 6.3.1.1494)
    
# For some reason, OSX likes it if I link these libraries before linking
# them all at the end.
target_link_libraries(qhullcpp qhull_p)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src)

set(library_SOURCES 
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPBoundaryUtils.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPConfigManager.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPConfigMessenger.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPCrystalGroup.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPDriftBoundaryProcess.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPDriftElectron.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPDriftHole.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPDriftRecombinationProcess.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPEigenSolver.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPElectrodeHit.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPElectrodeSensitivity.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPEnergyLimiter.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPEnergyPartition.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPEqEMField.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPFieldManager.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPInterpolator.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPInterValleyScattering.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPLocalElectroMagField.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPLukeScattering.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPMeshElectricField.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPMessengerBase.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPPhononBoundaryProcess.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPPhononKinTable.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPPhononKinematics.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPPhysics.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPPhysicsList.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPProcessUtils.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPSecondaryProduction.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPStackingAction.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPSurfaceProperty.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPTimeStepper.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPTrackInformation.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPTriLinearInterp.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPUnitsTable.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPUtils.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPVDriftProcess.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPVElectrodePattern.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4LatticeLogical.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4LatticeManager.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4LatticePhysical.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4LatticeReader.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4PhononDownconversion.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4PhononLong.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4PhononPolarization.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4PhononScattering.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4PhononTransFast.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4PhononTransSlow.cc
   ${CMAKE_CURRENT_SOURCE_DIR}/src/G4VPhononProcess.cc
   )

set(library_HEADERS
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPBlockData.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPBlockData.icc
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPBoundaryUtils.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPConfigManager.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPConfigMessenger.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPCrystalGroup.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPDriftBoundaryProcess.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPDriftElectron.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPDriftHole.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPDriftRecombinationProcess.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPEigenSolver.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPElectrodeHit.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPElectrodeSensitivity.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPEnergyLimiter.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPEnergyPartition.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPEqEMField.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPFieldManager.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPInterpolator.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPInterValleyScattering.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPLocalElectroMagField.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPLukeScattering.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPMatrix.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPMatrix.icc
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPMeshElectricField.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPMessengerBase.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPPhononBoundaryProcess.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPPhononKinTable.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPPhononKinematics.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPPhysics.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPPhysicsList.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPProcessSubType.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPProcessUtils.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPSecondaryProduction.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPStackingAction.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPSurfaceProperty.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPTimeStepper.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPTrackInformation.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPTriLinearInterp.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPUtils.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPVDriftProcess.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPVElectrodePattern.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4LatticeLogical.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4LatticeManager.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4LatticePhysical.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4LatticeReader.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4PhononDownconversion.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4PhononLong.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4PhononPolarization.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4PhononScattering.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4PhononTransFast.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4PhononTransSlow.hh
   ${CMAKE_CURRENT_SOURCE_DIR}/include/G4VPhononProcess.hh
   )

set(config_FILES
    ${PROJECT_SOURCE_DIR}/g4cmp_env.sh
    ${PROJECT_SOURCE_DIR}/g4cmp_env.csh
    ${PROJECT_SOURCE_DIR}/G4CMPOrdParamTable.txt
    )

add_library(G4cmp SHARED ${library_SOURCES})
target_link_libraries(G4cmp ${Geant4_LIBRARIES} qhullcpp)
set_target_properties(G4cmp PROPERTIES
                        COMPILE_DEFINITIONS "qh_QHpointer")
install(TARGETS G4cmp qhullcpp qhull_p DESTINATION lib)
install(FILES ${library_HEADERS} DESTINATION include/G4CMP)
install(FILES ${config_FILES} DESTINATION bin COMPONENT config)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/CrystalMaps DESTINATION CrystalMaps COMPONENT config)
