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
include_directories(${PROJECT_SOURCE_DIR}/library/qhull-2012.1/src)

set(charge_SOURCES 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ChargeActionInitialization.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ChargeDetectorConstruction.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ChargeElectrodeSensitivity.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ChargeFETDigitizerModule.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ChargeFETDigitizerMessenger.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ChargePrimaryGeneratorAction.cc
    )

set(fet_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ChargeFETDigitizerModule.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ChargeFETDigitizerMessenger.cc
    )

set(ionize_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ChargeActionInitialization.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ChargeDetectorConstruction.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ChargeElectrodeSensitivity.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ChargeFETDigitizerModule.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ChargeFETDigitizerMessenger.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ChargePrimaryGeneratorAction.cc
    )

set(charge_EPOT_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/EPotFiles/EPot_iZIP4_4V
    ${CMAKE_CURRENT_SOURCE_DIR}/EPotFiles/EPot_iZIP4_4V_small
    ${CMAKE_CURRENT_SOURCE_DIR}/EPotFiles/EPot_iZIP4_70V_small
    ${CMAKE_CURRENT_SOURCE_DIR}/EPotFiles/EPot_iZIP4_axissymmetric_shortQOS1
    )

set(charge_MACROS
    ${CMAKE_CURRENT_SOURCE_DIR}/run.mac
    ${CMAKE_CURRENT_SOURCE_DIR}/vis.mac
    ${CMAKE_CURRENT_SOURCE_DIR}/fast.mac
    ${CMAKE_CURRENT_SOURCE_DIR}/loop.mac
    ${CMAKE_CURRENT_SOURCE_DIR}/open.mac
    ${CMAKE_CURRENT_SOURCE_DIR}/play.mac
    ${CMAKE_CURRENT_SOURCE_DIR}/iZip4.mac
    ${CMAKE_CURRENT_SOURCE_DIR}/single.mac
    ${CMAKE_CURRENT_SOURCE_DIR}/openNew.mac
    ${CMAKE_CURRENT_SOURCE_DIR}/useRndm.mac
    ${CMAKE_CURRENT_SOURCE_DIR}/saveRndm.mac
    ${CMAKE_CURRENT_SOURCE_DIR}/testRefl.mac
    ${CMAKE_CURRENT_SOURCE_DIR}/shortQOS1.mac
    ${CMAKE_CURRENT_SOURCE_DIR}/testIonize.mac
    ${CMAKE_CURRENT_SOURCE_DIR}/drift_curve.mac
    ${CMAKE_CURRENT_SOURCE_DIR}/movie.loop
    )

set(fet_CONFIGS
    ${CMAKE_CURRENT_SOURCE_DIR}/config/ConstantsFET
    ${CMAKE_CURRENT_SOURCE_DIR}/config/FETTemplates
    ${CMAKE_CURRENT_SOURCE_DIR}/config/EpotRamoChan1
    ${CMAKE_CURRENT_SOURCE_DIR}/config/EpotRamoChan2
    ${CMAKE_CURRENT_SOURCE_DIR}/config/EpotRamoChan3
    ${CMAKE_CURRENT_SOURCE_DIR}/config/EpotRamoChan4
    )

add_executable(g4cmpCharge g4cmpCharge.cc ${charge_SOURCES})
target_link_libraries(g4cmpCharge G4cmp)
set_target_properties(g4cmpCharge PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib"
    INSTALL_RPATH_USE_LINK_RPATH TRUE
    BUILD_WITH_INSTALL_RPATH FALSE)

add_executable(g4cmpIonize g4cmpIonize.cc ${ionize_SOURCES})
target_link_libraries(g4cmpIonize G4cmp)
set_target_properties(g4cmpIonize PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib"
    INSTALL_RPATH_USE_LINK_RPATH TRUE
    BUILD_WITH_INSTALL_RPATH FALSE)

add_executable(g4cmpFETSim g4cmpFETSim.cc ${fet_SOURCES})
target_link_libraries(g4cmpFETSim G4cmp)
set_target_properties(g4cmpFETSim PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib"
    INSTALL_RPATH_USE_LINK_RPATH TRUE
    BUILD_WITH_INSTALL_RPATH FALSE)

install(TARGETS g4cmpCharge DESTINATION bin)
install(TARGETS g4cmpIonize DESTINATION bin)
install(TARGETS g4cmpFETSim DESTINATION bin)
install(FILES ${charge_EPOT_FILES} DESTINATION EPotFiles COMPONENT config)
install(FILES ${charge_MACROS} DESTINATION macros/g4cmpCharge)
install(FILES ${fet_CONFIGS} DESTINATION config/G4CMP/FETSim)
