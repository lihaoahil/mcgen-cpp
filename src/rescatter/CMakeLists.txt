add_library(mcgen_rescatter STATIC
  vertex1.cpp
  vertex2.cpp
  event.cpp
  writer_ascii.cpp
  writer_truth.cpp
)
target_link_libraries(mcgen_rescatter PUBLIC mcgen_core)

add_executable(mc_gen_rescatter main_gen_rescatter.cpp)
target_link_libraries(mc_gen_rescatter PRIVATE mcgen_rescatter)
