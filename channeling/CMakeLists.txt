 
include_directories(${PROJECT_SOURCE_DIR}/channeling/include)

set(channeling_SOURCES 
   ${PROJECT_SOURCE_DIR}/channeling/src/A01DriftChamber.cc
   ${PROJECT_SOURCE_DIR}/channeling/src/A01DriftChamberHit.cc
   ${PROJECT_SOURCE_DIR}/channeling/src/A01EventAction.cc
   ${PROJECT_SOURCE_DIR}/channeling/src/A01EventActionMessenger.cc
   ${PROJECT_SOURCE_DIR}/channeling/src/ChannelingParticleUserInfo.cc
   ${PROJECT_SOURCE_DIR}/channeling/src/DetectorConstruction.cc
   ${PROJECT_SOURCE_DIR}/channeling/src/DetectorConstructionMessenger.cc
   ${PROJECT_SOURCE_DIR}/channeling/src/G4LindhardPartition.cc
   ${PROJECT_SOURCE_DIR}/channeling/src/G4ScreenedNuclearRecoil.cc
   ${PROJECT_SOURCE_DIR}/channeling/src/PhysicsList.cc
   ${PROJECT_SOURCE_DIR}/channeling/src/PhysicsListMessenger.cc
   ${PROJECT_SOURCE_DIR}/channeling/src/PrimaryGeneratorAction.cc
   ${PROJECT_SOURCE_DIR}/channeling/src/PrimaryGeneratorActionMessenger.cc
   ${PROJECT_SOURCE_DIR}/channeling/src/ProcessChanneling.cc
   ${PROJECT_SOURCE_DIR}/channeling/src/ProcessChannelingSlow.cc
   ${PROJECT_SOURCE_DIR}/channeling/src/StackingAction.cc
   ${PROJECT_SOURCE_DIR}/channeling/src/StackingMessenger.cc
   ${PROJECT_SOURCE_DIR}/channeling/src/StepMax.cc
   ${PROJECT_SOURCE_DIR}/channeling/src/StepMaxMessenger.cc
   ${PROJECT_SOURCE_DIR}/channeling/src/TrackingAction.cc
   ${PROJECT_SOURCE_DIR}/channeling/src/XCrystalCharacteristicArray.cc
   ${PROJECT_SOURCE_DIR}/channeling/src/XCrystalIntegratedDensityHub.cc
   ${PROJECT_SOURCE_DIR}/channeling/src/XCrystalIntegratedDensityPlanar.cc
   ${PROJECT_SOURCE_DIR}/channeling/src/XCrystalPlanarMoliereElectricField.cc
   ${PROJECT_SOURCE_DIR}/channeling/src/XCrystalPlanarMoliereElectronDensity.cc
   ${PROJECT_SOURCE_DIR}/channeling/src/XCrystalPlanarMolierePotential.cc
   ${PROJECT_SOURCE_DIR}/channeling/src/XCrystalPlanarMoliereTempPotential.cc
   ${PROJECT_SOURCE_DIR}/channeling/src/XCrystalPlanarNucleiDensity.cc
   ${PROJECT_SOURCE_DIR}/channeling/src/XLatticManager3.cc
   ${PROJECT_SOURCE_DIR}/channeling/src/XLogicalAtomicLattice.cc
   ${PROJECT_SOURCE_DIR}/channeling/src/XLogicalAtomicLatticeDiamond.cc
   ${PROJECT_SOURCE_DIR}/channeling/src/XLogicalAtomicLatticeDiamondSingleAtom.cc
   ${PROJECT_SOURCE_DIR}/channeling/src/XLogicalBase.cc
   ${PROJECT_SOURCE_DIR}/channeling/src/XLogicalLattice.cc
   ${PROJECT_SOURCE_DIR}/channeling/src/XPhysicalLattice.cc
   ${PROJECT_SOURCE_DIR}/channeling/src/XUnitCell.cc
   ${PROJECT_SOURCE_DIR}/channeling/src/XVCrystalCharacteristic.cc
   ${PROJECT_SOURCE_DIR}/channeling/src/XVCrystalIntegratedDensity.cc
   ${PROJECT_SOURCE_DIR}/channeling/src/XVCrystalPlanarAnalytical.cc
   ${PROJECT_SOURCE_DIR}/channeling/src/XWrapperContinuousDiscreteProcess.cc
   ${PROJECT_SOURCE_DIR}/channeling/src/XWrapperDiscreteProcess.cc
   )
    add_executable(G4CMP_channeling g4cmpChannel.cc ${channeling_SOURCES})
    target_link_libraries(G4CMP_channeling ${Geant4_LIBRARIES} G4CMP_lib)
    set_target_properties(G4CMP_channeling PROPERTIES
                        COMPILE_DEFINITIONS "qh_QHpointer")
    set(CMAKE_CXX_FLAGS "-O2 -p -g")
    install(TARGETS G4CMP_channeling DESTINATION bin)