cmake_minimum_required(VERSION 3.1)

set(PROJECT_NAME cec)


project(${PROJECT_NAME} VERSION 1.0
             DESCRIPTION "Re-implementation of CECs benchmarks in C"
             LANGUAGES C
)

include(CTest)
enable_testing()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(EXECUTABLE ${PROJECT_NAME})
set(TEST ${PROJECT_NAME}-test)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
  ${PROJECT_SOURCE_SIR}/bin/${CMAKE_BUILD_TYPE})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

set(EXTRA_LIBS m)

find_path(CMOCKA_INCLUDE_DIR
  NAMES
  cmocka.h
  PATHS
  $(CMOCKA_INCLUDE_DIR)/include
)

find_library(CMOCKA_LIBRARY
  NAMES
  cmocka cmocka_shared
  PATHS
  ${CMOCKA_ROOT_DIR}/include
)

if (CMOCKA_LIBRARY)
  set(CMOCKA_LIBRARIES
    ${CMOCKA_LIBRARIES}
    ${CMOCKA_LIBRARY}
  )
endif (CMOCKA_LIBRARY)


include_directories(${CMOCKA_INCLUDE_DIR} include)

file(GLOB_RECURSE SRC_FILES "${PROJECT_SOURCE_DIR}/src/*.c")
file(GLOB_RECURSE INC_FILES "${PROJECT_SOURCE_DIR}/include/*.h")
file(GLOB_RECURSE TEST_FILES "${PROJECT_SOURCE_DIR}/test/*.c")

add_compile_options(-Wall -g -Wextra -pedantic -Werror)

add_executable(${EXECUTABLE} main.c ${SRC_FILES} ${INC_FILES})
target_link_libraries(${EXECUTABLE} ${CMOCKA_LIBRARIES} ${EXTRA_LIBS})
add_executable(${TEST} ${TEST_FILES} ${SRC_FILES} ${INC_FILES})
target_link_libraries(${TEST} ${CMOCKA_LIBRARIES} ${EXTRA_LIBS})

add_test(runTest "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${TEST}")


