add_library(mcgen_core STATIC
  fourvector.cpp
  rng.cpp
  defparser.cpp
  lineshape.cpp
  nbody.cpp
  beamgen.cpp
  dynamics.cpp
  targetgeom.cpp
  particleid.cpp
  generator.cpp
)
target_include_directories(mcgen_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
