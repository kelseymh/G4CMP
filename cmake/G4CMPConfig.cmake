# Find ourself. Very zen.
get_filename_component(_thisdir "${CMAKE_CURRENT_LIST_FILE}" PATH)

# Find all libraries included with G4CMP.
set(G4CMP_BASEDIR "${_thisdir}/..")
set(G4CMP_LIBDIR "${G4CMP_BASEDIR}/lib")

find_library(G4CMP_LIBRARY NAMES G4cmp HINTS ${G4CMP_LIBDIR})
find_library(QHCPP_LIBRARY NAMES qhullcpp HINTS ${G4CMP_LIBDIR})
find_library(QH_P_LIBRARY NAMES qhull_p HINTS ${G4CMP_LIBDIR})

# Set a bunch of useful variables.
set(G4CMP_INCLUDE_DIRS "${G4CMP_BASEDIR}/include/G4CMP")
set(G4CMP_LIBRARIES ${G4CMP_LIBRARY} ${QHCPP_LIBRARY} ${QH_P_LIBRARY})
set(G4CMP_CXX_FLAGS "-std=c++11")

# Users of G4CMP should add "include(${G4CMP_USE_FILE})" to CMakeLists.txt
set(G4CMP_USE_FILE ${G4CMP_BASEDIR}/cmake/UseG4CMP.cmake)
