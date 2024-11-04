cmake_minimum_required(VERSION 3.10)
project(som)
#add_definitions(-std=c++0x)
set(CMAKE_CXX_FLAGS='-std=c++17')
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_BUILD_TYPE Debug)

# Set the output directory for the build executables and libraries
#SET(som_RUNTIME_OUTPUT_DIRECTORY ${som_SOURCE_DIR}/bin CACHE PATH "Target for the binaries")
#SET(som_LIBRARY_OUTPUT_DIRECTORY ${som_SOURCE_DIR}/lib CACHE PATH "Target for the libraries")
#SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${som_LIBRARY_OUTPUT_DIRECTORY})
#SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${som_LIBRARY_OUTPUT_DIRECTORY})
#SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${som_RUNTIME_OUTPUT_DIRECTORY})
#SET(som_ENABLE_PCL 0)

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")

include_directories(include)
#include_directories(src)

#list(APPEND som_LINK_LIBS boost_iostreams boost_system boost_filesystem)

# Eigen library v. 3
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# rofl
find_package(rofl REQUIRED)
include_directories(${rofl_INCLUDE_DIR})

#if(${rofl_FOUND})
#  list(APPEND ${rofl_LIBRARIES})
#endif()

# Library files
#add_library(som
#  src/som2d.cpp
#)

add_library(roptlib_local
include/thirdparty/roptlib/manifolds_CFixedRankQ2F.cpp
include/thirdparty/roptlib/manifolds_CStiefel.cpp
include/thirdparty/roptlib/manifolds_CSymFixedRankQ.cpp
include/thirdparty/roptlib/manifolds_Element.cpp
include/thirdparty/roptlib/manifolds_Euclidean.cpp
include/thirdparty/roptlib/manifolds_FixedRankE.cpp
include/thirdparty/roptlib/manifolds_FixedRankQ2F.cpp
include/thirdparty/roptlib/manifolds_Grassmann.cpp
include/thirdparty/roptlib/manifolds_Manifold.cpp
include/thirdparty/roptlib/manifolds_MultiManifolds.cpp
include/thirdparty/roptlib/manifolds_Rotations.cpp #new
include/thirdparty/roptlib/manifolds_SmartSpace.cpp
include/thirdparty/roptlib/manifolds_SPDManifold.cpp
include/thirdparty/roptlib/manifolds_Sphere.cpp
# include/thirdparty/roptlib/manifolds_SphereTx.cpp
include/thirdparty/roptlib/manifolds_Stiefel.cpp
include/thirdparty/roptlib/manifolds_SymFixedRankQ.cpp
#others
include/thirdparty/roptlib/others_BlasLapackCppWrapper.cpp
include/thirdparty/roptlib/others_ForDebug.cpp
include/thirdparty/roptlib/others_MinPNormConHull.cpp
include/thirdparty/roptlib/others_randgen.cpp
include/thirdparty/roptlib/others_SparseMatrix.cpp
include/thirdparty/roptlib/others_Spline.cpp
# include/thirdparty/roptlib/others_Timer.cpp
include/thirdparty/roptlib/others_SparseBLAS_nist_spblas.cpp
#problems
include/thirdparty/roptlib/problems_CFRankQ2FBlindDecon2D.cpp
include/thirdparty/roptlib/problems_CSFRQPhaseRetrieval.cpp
include/thirdparty/roptlib/problems_CStieBrockett.cpp
include/thirdparty/roptlib/problems_EucQuadratic.cpp
include/thirdparty/roptlib/problems_FRankESparseApprox.cpp
include/thirdparty/roptlib/problems_FRankEWeightApprox.cpp
include/thirdparty/roptlib/problems_FRankQ2FMatCompletion.cpp
include/thirdparty/roptlib/problems_GrassRQ.cpp
include/thirdparty/roptlib/problems_juliaProblem.cpp
include/thirdparty/roptlib/problems_mexProblem.cpp
include/thirdparty/roptlib/problems_Problem.cpp
include/thirdparty/roptlib/problems_ProdStieSumBrockett.cpp
include/thirdparty/roptlib/problems_SampleSom.cpp #som
include/thirdparty/roptlib/problems_SFRQLyapunov.cpp
include/thirdparty/roptlib/problems_SPDKarcherMean.cpp
include/thirdparty/roptlib/problems_SphereConvexHull.cpp
include/thirdparty/roptlib/problems_SphereSparsestVector.cpp
include/thirdparty/roptlib/problems_SphereTxRQ.cpp
include/thirdparty/roptlib/problems_StieBrockett.cpp
include/thirdparty/roptlib/problems_StieSPCA.cpp
include/thirdparty/roptlib/manifolds_SphereTx.cpp
#solvers
include/thirdparty/roptlib/solvers_AManPG.cpp
include/thirdparty/roptlib/solvers_LRBFGS.cpp
include/thirdparty/roptlib/solvers_LRBFGSSub.cpp
include/thirdparty/roptlib/solvers_LRTRSR1.cpp
include/thirdparty/roptlib/solvers_ManPG.cpp
include/thirdparty/roptlib/solvers_RBFGS.cpp
include/thirdparty/roptlib/solvers_RBFGSSub.cpp
include/thirdparty/roptlib/solvers_RBroydenFamily.cpp
include/thirdparty/roptlib/solvers_RCG.cpp
include/thirdparty/roptlib/solvers_RGS.cpp
include/thirdparty/roptlib/solvers_RNewton.cpp
include/thirdparty/roptlib/solvers_RSD.cpp
include/thirdparty/roptlib/solvers_RTRNewton.cpp
include/thirdparty/roptlib/solvers_RTRSD.cpp
include/thirdparty/roptlib/solvers_RTRSR1.cpp
include/thirdparty/roptlib/solvers_RWRBFGS.cpp
include/thirdparty/roptlib/solvers_Solvers.cpp
include/thirdparty/roptlib/solvers_SolversNSM.cpp
include/thirdparty/roptlib/solvers_SolversNSMPGLS.cpp
include/thirdparty/roptlib/solvers_SolversNSMSub.cpp
include/thirdparty/roptlib/solvers_SolversNSMSubLS.cpp
include/thirdparty/roptlib/solvers_SolversSM.cpp
include/thirdparty/roptlib/solvers_SolversSMLS.cpp
include/thirdparty/roptlib/solvers_SolversSMTR.cpp
include/thirdparty/roptlib/others_Timer.cpp
)

add_executable(test_roptlib_som test_roptlib_som.cpp)
target_link_libraries(test_roptlib_som ${rofl_LIBRARIES} roptlib_local m blas lapack boost_iostreams boost_system boost_filesystem)

# Option "make install": copy binaries
#INSTALL(TARGETS som
#  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
#  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
#  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
#)

# Option "make install": copy headers
#FILE(GLOB headers "${CMAKE_CURRENT_SOURCE_DIR}/include/som/*.h")
#FILE(GLOB headers_3rdparty "${CMAKE_CURRENT_SOURCE_DIR}/include/som/thirdparty/*.h")
#INSTALL(FILES ${headers} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/som/)
#INSTALL(FILES ${headers_3rdparty} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/som/thirdparty/)

# Option "make install": copy cmake script
#FILE(GLOB cmake_script "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/somConfig.cmake")
#message(STATUS "cmake_script " ${cmake_script})
#INSTALL(FILES ${cmake_script} DESTINATION ${CMAKE_INSTALL_PREFIX}/share/som/)

# uninstall target
#if(NOT TARGET uninstall)
#  configure_file(
#    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/cmake_uninstall.cmake.in"
#    "${CMAKE_CURRENT_BINARY_DIR}/cmake_modules/cmake_uninstall.cmake"
#    IMMEDIATE @ONLY)

#  add_custom_target(uninstall
#    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_modules/cmake_uninstall.cmake)
#endif()
