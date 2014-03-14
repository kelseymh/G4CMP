include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
set(charge_SOURCES 
    ${PROJECT_SOURCE_DIR}/charge/src/CDMS_iZip4_Field.cc
    ${PROJECT_SOURCE_DIR}/charge/src/ChargeDetectorConstruction.cc
    ${PROJECT_SOURCE_DIR}/charge/src/ChargeEMField.cc
    ${PROJECT_SOURCE_DIR}/charge/src/ChargePrimaryGeneratorAction.cc
    ${PROJECT_SOURCE_DIR}/charge/src/EqEMFieldXtal.cc
    ${PROJECT_SOURCE_DIR}/charge/src/FET.cc
    ${PROJECT_SOURCE_DIR}/charge/src/FETMessenger.cc
    ${PROJECT_SOURCE_DIR}/charge/src/MatWriter.cc
    ${PROJECT_SOURCE_DIR}/charge/src/XtalFieldManager.cc
    )
add_executable(G4CMP_charge g4cmpCharge.cc ${charge_SOURCES})
target_link_libraries(G4CMP_charge ${Geant4_LIBRARIES} G4CMP_lib)
set_target_properties(G4CMP_charge PROPERTIES
                        COMPILE_DEFINITIONS "qh_QHpointer")
set(CMAKE_CXX_FLAGS "-O2 -p -g")
install(TARGETS G4CMP_charge DESTINATION bin)
