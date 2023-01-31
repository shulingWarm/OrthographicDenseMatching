file(GLOB TESTS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

# A semi-colon separated list of test sources that should not be automatically built with boost unit test
set(SPECIAL_TESTS "portability_test.cpp")

# Build the portability test only if we are on a 64-bit machine (void* is 8 bytes)
if((${CMAKE_SIZEOF_VOID_P} EQUAL 8) AND (NOT SKIP_PORTABILITY_TEST))

  add_test(NAME portability_test COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/run_portability_test.sh")

endif()

# Build all of the non-special tests and link against the boost unit test framework
foreach(TEST_SOURCE ${TESTS})

  string(REPLACE ".cpp" "" TEST_TARGET "${TEST_SOURCE}")
  set(TEST_TARGET "test_${TEST_TARGET}")

  # Check to see if our target is listed in "SPECIAL_TESTS"
  list(FIND SPECIAL_TESTS "${TEST_SOURCE}" IS_SPECIAL_TEST)

endforeach()

# Add the valgrind target
add_custom_target(valgrind
  COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/run_valgrind.sh")

# Add the coverage target
add_custom_target(coverage)
add_custom_command(TARGET coverage
  COMMAND "${CMAKE_SOURCE_DIR}/scripts/updatecoverage.sh" ${CMAKE_SOURCE_DIR}
  WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/coverage")

# add tests to coverage
foreach(TEST_SOURCE ${TESTS})
  string(REPLACE ".cpp" "" COVERAGE_TARGET "${TEST_SOURCE}")
  set(COVERAGE_TARGET "coverage_${COVERAGE_TARGET}")

  # Check to see if our target is listed in "SPECIAL_TESTS"
  list(FIND SPECIAL_TESTS "${TEST_SOURCE}" IS_SPECIAL_TEST)

  if(IS_SPECIAL_TEST EQUAL -1)
    add_dependencies(coverage ${COVERAGE_TARGET})
  endif()
endforeach()

if(NOT CMAKE_VERSION VERSION_LESS 3.0)
  add_test(test_cmake_config_module ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake-config-module.cmake)
endif()
