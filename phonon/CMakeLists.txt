include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

set(phonon_SOURCES 
   ${PROJECT_SOURCE_DIR}/phonon/src/PhononDetectorConstruction.cc
   ${PROJECT_SOURCE_DIR}/phonon/src/PhononPrimaryGeneratorAction.cc
   )
    
add_executable(G4CMP_phonon g4cmpPhonon.cc ${phonon_SOURCES})
target_link_libraries(G4CMP_phonon ${Geant4_LIBRARIES} G4CMP_lib)
set_target_properties(G4CMP_phonon PROPERTIES
                        COMPILE_DEFINITIONS "qh_QHpointer")
set(CMAKE_CXX_FLAGS "-O2 -p -g")
install(TARGETS G4CMP_phonon DESTINATION bin)