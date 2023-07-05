#----------------------------------------------------------------------------
# G4CMP-specific command-line options (i.e., -DXYZ=ON or -DXYZ=OFF
#

set(G4CMP_DEBUG "" CACHE STRING "Enable debugging blocks and verbosity")
set_property(CACHE G4CMP_DEBUG PROPERTY STRINGS "" 0 1 2 3 4)

option(G4CMPTLI_DEBUG "Enable debugging of TriLinearInterp" OFF)

#----------------------------------------------------------------------------
# Sanitize Multithreaded code
# Use -DG4CMP_USE_SANITIZER=ON to enable "-fsanitize=thread" on compilation.
# Use -DG4CMP_SANITIZER_TYPE=thread (or memory, address, leak)
#
# Adapted rom Jonathan Madsen, LBL
#
option(G4CMP_USE_SANITIZER "Enable -fsanitize=\${G4CMP_SANITIZER_TYPE}" OFF)

# the type of sanitizer to use, e.g. thread, memory, address, leak
set(G4CMP_SANITIZER_TYPE "thread" CACHE STRING "-fsanitize=<TYPE>")
set_property(CACHE G4CMP_SANITIZER_TYPE PROPERTY
	     STRINGS "" "thread" "memory" "address" "leak")

if(G4CMP_USE_SANITIZER)
   set(G4CMP_FLAGS_SANITIZE "-fno-omit-frame-pointer" "-fsanitize=${G4CMP_SANITIZER_TYPE}")
endif()

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
# Qhull Compilation stuff -- sources are under current (library/) directory
#
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhull)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1/src/libqhullcpp)

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
)

if(USE_GEANT4_STATIC_LIBS)
    add_library(qhull_p STATIC ${libqhull_SOURCES})
else()
    add_library(qhull_p SHARED ${libqhull_SOURCES})
    if(APPLE)
        set_property(TARGET qhull_p
	    APPEND PROPERTY INSTALL_RPATH "@loader_path")
    endif()
endif()
set_target_properties(qhull_p PROPERTIES
    COMPILE_DEFINITIONS "qh_QHpointer"
    VERSION 6.3.1.1494)

if(USE_GEANT4_STATIC_LIBS)
    add_library(qhullcpp STATIC ${libqhullcpp_SOURCES})
else()
    add_library(qhullcpp SHARED ${libqhullcpp_SOURCES})
    if(APPLE)
        set_property(TARGET qhullcpp
	    APPEND PROPERTY INSTALL_RPATH "@loader_path")
    endif()
endif()
set_target_properties(qhullcpp PROPERTIES
    COMPILE_DEFINITIONS "qh_QHpointer"
    VERSION 6.3.1.1494)

add_library(qhullstatic_p STATIC ${libqhull_SOURCES})
set_target_properties(qhullstatic_p PROPERTIES
    COMPILE_DEFINITIONS "qh_QHpointer"
    VERSION 6.3.1.1494)
    
# For some reason, OSX likes it if I link these libraries before linking
# them all at the end.
target_link_libraries(qhullcpp qhull_p qhullstatic_p)

add_definitions(-DG4DIGI_ALLOC_EXPORT)		## Needed for hits collection

set(library_SOURCES 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPAnharmonicDecay.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPBiLinearInterp.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPBoundaryUtils.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPChargeCloud.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPConfigManager.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPConfigMessenger.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPCrystalGroup.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPDownconversionRate.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPDriftBoundaryProcess.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPDriftElectron.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPDriftHole.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPDriftTrapIonization.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPDriftRecombinationProcess.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPDriftTrackInfo.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPDriftTrappingProcess.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPEigenSolver.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPElectrodeHit.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPElectrodeSensitivity.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPEnergyPartition.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPEqEMField.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPFanoBinomial.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPFieldManager.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPFieldUtils.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPGeometryUtils.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPGlobalLocalTransformStore.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPHitMerging.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPIVRateLinear.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPIVRateQuadratic.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPInterValleyRate.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPInterValleyScattering.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPInterpolator.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPKaplanQP.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPLewinSmithNIEL.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPLindhardNIEL.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPLocalElectroMagField.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPLogicalBorderSurface.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPLogicalSkinSurface.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPLukeEmissionRate.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPLukeScattering.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPMeshElectricField.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPPartitionData.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPPartitionSummary.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPPhononBoundaryProcess.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPPhononElectrode.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPPhononKinTable.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPPhononKinematics.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPPhononScatteringRate.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPPhononTrackInfo.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPPhysics.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPPhysicsList.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPProcessUtils.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPSecondaryProduction.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPSecondaryUtils.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPStackingAction.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPStepAccumulator.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPSurfaceProperty.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPTimeStepper.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPTrackLimiter.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPTrackUtils.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPTriLinearInterp.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPUnitsTable.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPUtils.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPVDriftProcess.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPVElectrodePattern.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPVMeshInterpolator.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPVProcess.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4CMPVTrackInfo.cc
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
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4VNIELPartition.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/G4VPhononProcess.cc
    )
 
set(library_HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPAnharmonicDecay.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPBiLinearInterp.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPBlockData.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPBlockData.icc
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPBoundaryUtils.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPChargeCloud.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPConfigManager.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPConfigMessenger.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPCrystalGroup.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPDownconversionRate.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPDriftBoundaryProcess.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPDriftElectron.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPDriftHole.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPDriftTrapIonization.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPDriftRecombinationProcess.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPDriftTrackInfo.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPDriftTrappingProcess.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPEigenSolver.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPElectrodeHit.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPElectrodeSensitivity.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPEnergyPartition.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPEqEMField.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPFanoBinomial.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPFanoBinomial.icc
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPFieldManager.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPFieldUtils.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPGeometryUtils.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPGlobalLocalTransformStore.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPHitMerging.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPIVRateLinear.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPIVRateQuadratic.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPInterValleyRate.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPInterValleyScattering.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPInterpolator.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPKaplanQP.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPLewinSmithNIEL.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPLindhardNIEL.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPLocalElectroMagField.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPLogicalBorderSurface.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPLogicalSkinSurface.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPLukeEmissionRate.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPLukeScattering.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPMatrix.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPMatrix.icc
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPMeshElectricField.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPPartitionData.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPPartitionSummary.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPPhononBoundaryProcess.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPPhononElectrode.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPPhononKinTable.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPPhononKinematics.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPPhononScatteringRate.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPPhononTrackInfo.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPPhysics.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPPhysicsList.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPProcessSubType.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPProcessUtils.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPSecondaryProduction.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPSecondaryUtils.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPStackingAction.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPStepAccumulator.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPSurfaceProperty.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPTimeStepper.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPTrackLimiter.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPTrackUtils.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPTrackUtils.icc
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPTriLinearInterp.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPUnitsTable.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPUtils.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPVDriftProcess.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPVElectrodePattern.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPVMeshInterpolator.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPVProcess.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPVScatteringRate.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4CMPVTrackInfo.hh
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
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4VNIELPartition.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4VPhononProcess.hh
    )

set(config_FILES
    ${PROJECT_SOURCE_DIR}/g4cmp_env.sh
    ${PROJECT_SOURCE_DIR}/g4cmp_env.csh
    ${PROJECT_SOURCE_DIR}/g4cmp.gmk
    ${PROJECT_SOURCE_DIR}/FindGeant4.cmake
    )

if(USE_GEANT4_STATIC_LIBS)
    add_library(G4cmp STATIC ${library_SOURCES})
else()
    add_library(G4cmp SHARED ${library_SOURCES})
    if(APPLE)
        set_property(TARGET G4cmp
	    APPEND PROPERTY INSTALL_RPATH "@loader_path")
    endif()
endif()

target_link_libraries(G4cmp PUBLIC ${Geant4_LIBRARIES} qhullcpp)

set(LibDefs "qh_QHpointer")
if(NOT G4CMP_DEBUG STREQUAL "")
    set(LibDefs "${LibDefs};G4CMP_DEBUG=${G4CMP_DEBUG}")
endif()
if(G4CMPTLI_DEBUG)
    set(LibDefs "${LibDefs};G4CMPTLI_DEBUG=1")
endif()
if(Geant4_builtin_clhep_FOUND)
    SET(LibDefs "${LibDefs};G4LIB_USE_CLHEP=1")
endif()

target_compile_definitions(G4cmp PUBLIC ${LibDefs})

if(G4CMP_USE_SANITIZER)
    target_compile_options(G4cmp INTERFACE ${G4CMP_FLAGS_SANITIZE})
    set_property(TARGET G4cmp PROPERTY
        INTERFACE_LINK_OPTIONS ${G4CMP_FLAGS_SANITIZE})
endif()

install(TARGETS G4cmp qhullcpp qhull_p qhullstatic_p DESTINATION lib)
install(FILES ${library_HEADERS} DESTINATION include/G4CMP)
install(FILES ${config_FILES} DESTINATION share/G4CMP COMPONENT config)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/CrystalMaps DESTINATION share/G4CMP COMPONENT config)
