

cmake_minimum_required(VERSION 3.10)
project(qr_unique_sizeless)
#add_definitions(-std=c++0x)
set(CMAKE_CXX_FLAGS='-std=c++17')
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_BUILD_TYPE Debug)

include_directories(include_matlab)
#include_directories(src)

#list(APPEND som_LINK_LIBS boost_iostreams boost_system boost_filesystem)

# Eigen library v. 3
find_package(Eigen3 3.4 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

add_library(qru
qr_unique.cpp
qr_unique_initialize.cpp
qr_unique_terminate.cpp
qr.cpp
diag.cpp
#sign.cpp
#xgemv.cpp
#xgerc.cpp
xnrm2.cpp
#rtGetInf.cpp
#rtGetNaN.cpp
rt_nonfinite.cpp
)

add_executable(main main.cpp)
target_link_libraries(main qru)