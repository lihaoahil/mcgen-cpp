# Unit tests are added per-module as they are ported.
add_executable(test_defparser test_defparser.cpp)
target_link_libraries(test_defparser PRIVATE mcgen_core)

add_executable(test_kinematics test_kinematics.cpp)
target_link_libraries(test_kinematics PRIVATE mcgen_core)
add_test(NAME kinematics COMMAND test_kinematics)
