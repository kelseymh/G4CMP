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
    ${PROJECT_SOURCE_DIR}/charge/src/ChargeEMField.cc
    ${PROJECT_SOURCE_DIR}/charge/src/ChargeFieldMessenger.cc
    ${PROJECT_SOURCE_DIR}/charge/src/ChargePrimaryGeneratorAction.cc
    ${PROJECT_SOURCE_DIR}/charge/src/FET.cc
    ${PROJECT_SOURCE_DIR}/charge/src/FETMessenger.cc
    ${PROJECT_SOURCE_DIR}/charge/src/MatWriter.cc
    )
add_executable(g4cmpCharge g4cmpCharge.cc ${charge_SOURCES})
target_link_libraries(g4cmpCharge G4cmp ${Geant4_LIBRARIES} qhullcpp qhullstatic_p)
set_target_properties(g4cmpCharge PROPERTIES
                        COMPILE_DEFINITIONS "qh_QHpointer")
set(CMAKE_CXX_FLAGS "-O2 -p -g")
install(TARGETS g4cmpCharge DESTINATION bin)
