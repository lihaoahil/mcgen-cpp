add_library(mcgen_analysis STATIC
  rootanalyzer.cpp
)
target_include_directories(mcgen_analysis PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_link_libraries(mcgen_analysis
  PUBLIC mcgen_core
  PRIVATE ROOT::Hist ROOT::RIO)

# Hao's ROOT 6.32 build was compiled with C++17; ROOT headers are only
# included inside this target, so pin its standard to match the ROOT build.
set_target_properties(mcgen_analysis PROPERTIES
  CXX_STANDARD 17
  CXX_STANDARD_REQUIRED ON)
