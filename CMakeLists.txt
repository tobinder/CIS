project(CIS)
cmake_minimum_required(VERSION 2.6)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)

find_package(Vigra   REQUIRED)
find_package(X11     REQUIRED)
find_package(Threads REQUIRED)
find_package(Boost   REQUIRED)
find_package(HDF5    REQUIRED)
find_package(OpenMP  REQUIRED)

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${VIGRA_INCLUDE_DIR}
    ${BOOST_INCLUDE_DIR}
    ${HDF5_INCLUDE_DIR}
)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
add_executable(cis main.cpp)
target_link_libraries(cis ${VIGRA_IMPEX_LIBRARY} ${X11_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${HDF5_LIBRARY} ${HDF5_HL_LIBRARY})
