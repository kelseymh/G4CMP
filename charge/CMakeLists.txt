#----------------------------------------------------------------------------
# Need to include G4CMP directories before G4.
#
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/library/include)

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
# Grab Qhull files needed to link against it.
# 
include_directories(${PROJECT_SOURCE_DIR}/qhull-2012.1/src)

set(charge_SOURCES 
    ${PROJECT_SOURCE_DIR}/charge/src/ChargeActionInitialization.cc
    ${PROJECT_SOURCE_DIR}/charge/src/ChargeDetectorConstruction.cc
    ${PROJECT_SOURCE_DIR}/charge/src/ChargeElectrodeSensitivity.cc
    ${PROJECT_SOURCE_DIR}/charge/src/ChargeFETDigitizerModule.cc
    ${PROJECT_SOURCE_DIR}/charge/src/ChargeFETDigitizerMessenger.cc
    ${PROJECT_SOURCE_DIR}/charge/src/ChargePrimaryGeneratorAction.cc
    )

set(fet_SOURCES
    ${PROJECT_SOURCE_DIR}/charge/src/ChargeFETDigitizerModule.cc
    ${PROJECT_SOURCE_DIR}/charge/src/ChargeFETDigitizerMessenger.cc
    )


set(charge_EPOT_FILES
    ${PROJECT_SOURCE_DIR}/charge/EPotFiles/EPot_iZIP4_axissymmetric_shortQOS1
    ${PROJECT_SOURCE_DIR}/charge/EPotFiles/Epot_iZip4_small_shortQOS1
    ${PROJECT_SOURCE_DIR}/charge/EPotFiles/Epot_iZip4
    ${PROJECT_SOURCE_DIR}/charge/EPotFiles/Epot_iZip4_small
    ${PROJECT_SOURCE_DIR}/charge/EPotFiles/EPot_iZip4_70V_COMSOL5
    )

set(charge_MACROS
    ${PROJECT_SOURCE_DIR}/charge/run.mac
    ${PROJECT_SOURCE_DIR}/charge/vis.mac
    ${PROJECT_SOURCE_DIR}/charge/fast.mac
    ${PROJECT_SOURCE_DIR}/charge/loop.mac
    ${PROJECT_SOURCE_DIR}/charge/open.mac
    ${PROJECT_SOURCE_DIR}/charge/play.mac
    ${PROJECT_SOURCE_DIR}/charge/iZip4.mac
    ${PROJECT_SOURCE_DIR}/charge/single.mac
    ${PROJECT_SOURCE_DIR}/charge/openNew.mac
    ${PROJECT_SOURCE_DIR}/charge/useRndm.mac
    ${PROJECT_SOURCE_DIR}/charge/saveRndm.mac
    ${PROJECT_SOURCE_DIR}/charge/testRefl.mac
    ${PROJECT_SOURCE_DIR}/charge/shortQOS1.mac
    ${PROJECT_SOURCE_DIR}/charge/testIonize.mac
    ${PROJECT_SOURCE_DIR}/charge/drift_curve.mac
    ${PROJECT_SOURCE_DIR}/charge/movie.loop
    )

set(fet_CONFIGS
    ${PROJECT_SOURCE_DIR}/charge/config/ConstantsFET
    ${PROJECT_SOURCE_DIR}/charge/config/FETTemplates
    ${PROJECT_SOURCE_DIR}/charge/config/EpotRamoChan1
    ${PROJECT_SOURCE_DIR}/charge/config/EpotRamoChan2
    ${PROJECT_SOURCE_DIR}/charge/config/EpotRamoChan3
    ${PROJECT_SOURCE_DIR}/charge/config/EpotRamoChan4
    )

add_executable(g4cmpCharge g4cmpCharge.cc ${charge_SOURCES})
target_link_libraries(g4cmpCharge G4cmp ${Geant4_LIBRARIES} qhullcpp qhullstatic_p)

add_executable(g4cmpFETSim g4cmpFETSim.cc ${fet_SOURCES})
target_link_libraries(g4cmpFETSim G4cmp ${Geant4_LIBRARIES} qhullcpp qhullstatic_p)

install(TARGETS g4cmpCharge DESTINATION bin)
install(TARGETS g4cmpFETSim DESTINATION bin)
install(FILES ${charge_EPOT_FILES} DESTINATION EPotFiles COMPONENT config)
install(FILES ${charge_MACROS} DESTINATION macros/g4cmpCharge)
install(FILES ${fet_CONFIGS} DESTINATION config/G4CMP/FETSim)
