cmake_minimum_required(VERSION 3.16)
project(FreeFlow)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -m64 -O3")

find_package(Boost 1.7 REQUIRED)

include_directories(${GMSH_DIR}/include ${Boost_INCLUDE_DIR} ${EIGEN_DIR})
link_directories(${GMSH_DIR}/lib)

add_executable(FreeFlow
        src/FreeFlow.cpp
        src/include/EigenType.h
        src/Mesh/MeshHandler.cpp src/Mesh/MeshHandler.h
        src/include/CoordinateSystems.h
        src/include/Data.h
        src/Solver/FiniteVolume.cpp src/Solver/FiniteVolume.h
        src/System/NS.cpp src/System/NS.h src/System/VarName.h
        src/Schemes/SSPRK.cpp src/Schemes/SSPRK.h
        src/Schemes/Interpolator.cpp src/Schemes/Interpolator.h src/Solver/FiniteVolume.cpp src/Solver/FiniteVolume.h src/System/k_Omega.cpp src/System/k_Omega.h)

target_link_libraries(FreeFlow gmsh ${Boost_LIBRARIES})