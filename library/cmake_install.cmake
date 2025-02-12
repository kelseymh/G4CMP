# Install script for directory: /sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP-457_install")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libG4cmp.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libG4cmp.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libG4cmp.so"
         RPATH "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP-457_install/lib:/usr/local/lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/libG4cmp.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libG4cmp.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libG4cmp.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libG4cmp.so"
         OLD_RPATH "/usr/local/lib:/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library::::::::"
         NEW_RPATH "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP-457_install/lib:/usr/local/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libG4cmp.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libqhullcpp.so.6.3.1.1494" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libqhullcpp.so.6.3.1.1494")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libqhullcpp.so.6.3.1.1494"
         RPATH "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP-457_install/lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/libqhullcpp.so.6.3.1.1494")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libqhullcpp.so.6.3.1.1494" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libqhullcpp.so.6.3.1.1494")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libqhullcpp.so.6.3.1.1494"
         OLD_RPATH "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library::::::::"
         NEW_RPATH "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP-457_install/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libqhullcpp.so.6.3.1.1494")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libqhullcpp.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libqhullcpp.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libqhullcpp.so"
         RPATH "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP-457_install/lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/libqhullcpp.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libqhullcpp.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libqhullcpp.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libqhullcpp.so"
         OLD_RPATH "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library::::::::"
         NEW_RPATH "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP-457_install/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libqhullcpp.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libqhull_p.so.6.3.1.1494" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libqhull_p.so.6.3.1.1494")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libqhull_p.so.6.3.1.1494"
         RPATH "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP-457_install/lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/libqhull_p.so.6.3.1.1494")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libqhull_p.so.6.3.1.1494" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libqhull_p.so.6.3.1.1494")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libqhull_p.so.6.3.1.1494"
         OLD_RPATH "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
         NEW_RPATH "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP-457_install/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libqhull_p.so.6.3.1.1494")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libqhull_p.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libqhull_p.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libqhull_p.so"
         RPATH "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP-457_install/lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/libqhull_p.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libqhull_p.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libqhull_p.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libqhull_p.so"
         OLD_RPATH "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
         NEW_RPATH "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP-457_install/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libqhull_p.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/libqhullstatic_p.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/G4CMP" TYPE FILE FILES
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPAnharmonicDecay.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPBiLinearInterp.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPBlockData.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPBlockData.icc"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPBoundaryUtils.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPChargeCloud.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPConfigManager.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPConfigMessenger.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPCrystalGroup.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPDownconversionRate.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPDriftBoundaryProcess.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPDriftElectron.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPDriftHole.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPDriftTrapIonization.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPDriftRecombinationProcess.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPDriftTrackInfo.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPDriftTrappingProcess.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPEigenSolver.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPElectrodeHit.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPElectrodeSensitivity.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPEnergyPartition.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPEqEMField.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPFanoBinomial.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPFanoBinomial.icc"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPFieldManager.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPFieldUtils.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPGeometryUtils.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPGlobalLocalTransformStore.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPHitMerging.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPIVRateLinear.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPIVRateQuadratic.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPImpactTunlNIEL.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPInterValleyRate.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPInterValleyScattering.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPInterpolator.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPKaplanQP.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPLewinSmithNIEL.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPLindhardNIEL.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPLocalElectroMagField.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPLogicalBorderSurface.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPLogicalSkinSurface.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPLukeEmissionRate.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPLukeScattering.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPMatrix.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPMatrix.icc"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPMeshElectricField.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPPartitionData.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPPartitionSummary.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPPhononBoundaryProcess.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPPhononElectrode.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPPhononKinTable.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPPhononKinematics.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPPhononScatteringRate.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPPhononTrackInfo.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPPhysics.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPPhysicsList.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPProcessSubType.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPProcessUtils.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPSarkisNIEL.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPSecondaryProduction.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPSecondaryUtils.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPStackingAction.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPStepAccumulator.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPSurfaceProperty.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPTimeStepper.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPTrackLimiter.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPTrackUtils.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPTrackUtils.icc"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPTriLinearInterp.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPUnitsTable.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPUtils.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPVDriftProcess.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPVElectrodePattern.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPVMeshInterpolator.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPVProcess.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPVScatteringRate.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4CMPVTrackInfo.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4LatticeLogical.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4LatticeManager.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4LatticePhysical.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4LatticeReader.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4PhononDownconversion.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4PhononLong.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4PhononPolarization.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4PhononScattering.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4PhononTransFast.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4PhononTransSlow.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4VNIELPartition.hh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/library/include/G4VPhononProcess.hh"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xconfigx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/G4CMP" TYPE FILE FILES
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/g4cmp_env.sh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/g4cmp_env.csh"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/g4cmp.gmk"
    "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/FindGeant4.cmake"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xconfigx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/G4CMP" TYPE DIRECTORY FILES "/sdf/home/d/dhenis/soft/packages/G4CMP/G4CMP/CrystalMaps")
endif()

