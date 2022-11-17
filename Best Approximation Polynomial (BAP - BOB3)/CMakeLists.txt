cmake_minimum_required (VERSION 3.8 FATAL_ERROR)

project ("Best Approximation Polynomial (BAP - BOB3)")
find_package(Eigen3 REQUIRED)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bld)

add_executable (${PROJECT_NAME} "Best Approximation Polynomial (BAP - BOB3).cpp")

target_include_directories(${PROJECT_NAME} PRIVATE ${eigen3_dir})

target_link_libraries(${PROJECT_NAME} PRIVATE ${eigen3_library_dir})