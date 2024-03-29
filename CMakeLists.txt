cmake_minimum_required(VERSION 3.5)
project(MOproject)

enable_language(CXX)
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_FLAGS "-fopenmp")

find_package(HDF5 REQUIRED COMPONENTS C CXX)
find_package(Boost REQUIRED COMPONENTS filesystem)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

include_directories(${PROJECT_SOURCE_DIR}/include ${HDF5_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR} ${Boost_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)
add_executable(MO main.cpp ${sources} ${headers})
target_link_libraries(MO ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES} ${Boost_LIBRARIES})

install(TARGETS MO DESTINATION bin)
