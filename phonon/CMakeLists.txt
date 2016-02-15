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

set(phonon_SOURCES 
   ${PROJECT_SOURCE_DIR}/phonon/src/PhononActionInitialization.cc
   ${PROJECT_SOURCE_DIR}/phonon/src/PhononDetectorConstruction.cc
   ${PROJECT_SOURCE_DIR}/phonon/src/PhononPrimaryGeneratorAction.cc
   )
    
add_executable(g4cmpPhonon g4cmpPhonon.cc ${phonon_SOURCES})
target_link_libraries(g4cmpPhonon G4cmp ${Geant4_LIBRARIES} qhullcpp qhullstatic_p)
set(CMAKE_CXX_FLAGS "-std=c++11 -O2 -Wall")
install(TARGETS g4cmpPhonon DESTINATION bin)
