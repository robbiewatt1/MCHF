cmake_minimum_required(VERSION 3.5)
project(MOproject)

enable_language(CXX)
set (CMAKE_CXX_STANDARD 11)

#Flags
SET(GCC_COVERAGE_COMPILE_FLAGS "-ffast-math")
add_definitions(${GCC_COVERAGE_COMPILE_FLAGS})

find_package(HDF5 REQUIRED COMPONENTS C CXX)

include_directories(${PROJECT_SOURCE_DIR}/include ${HDF5_INCLUDE_DIRS})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(MO main.cpp ${sources} ${headers})
target_link_libraries(MO ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES})
install(TARGETS MO DESTINATION bin)
