cmake_minimum_required(VERSION 3.20)
project(mc_gen_cpp CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()

add_compile_options(-Wall -Wextra)

# Core generator library: pure C++20, no external dependencies.
add_subdirectory(src/core)

# Analysis library: requires ROOT. ROOT 6.32 at /Users/haoli/CERN_ROOT/root_build
# was built with C++17; keep all ROOT includes confined to src/analysis.
option(MCGEN_WITH_ROOT "Build the ROOT-based analysis/histogramming module" ON)
if(MCGEN_WITH_ROOT)
  find_package(ROOT QUIET HINTS $ENV{ROOTSYS} /Users/haoli/CERN_ROOT/root_build)
  if(ROOT_FOUND)
    message(STATUS "ROOT ${ROOT_VERSION} found at ${ROOT_DIR}")
    add_subdirectory(src/analysis)
  else()
    message(WARNING "ROOT not found - building generator without analysis module")
    set(MCGEN_WITH_ROOT OFF)
  endif()
endif()

add_subdirectory(src/app)

enable_testing()
add_subdirectory(tests)
