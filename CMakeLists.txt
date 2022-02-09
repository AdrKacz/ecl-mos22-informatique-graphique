# Minimum version of CMake required to build this project
cmake_minimum_required(VERSION 3.0)

# Set specific std
add_definitions(-std=c++20)

# How to set option -Ipthread for g++ ??

# Name of the project
project(InformatiqueGraphique)

# Create a few variables for the folder names, so they are easier to rename in
# the future
set(SRC_DIR src)

# Add all the source files needed to build the executable
add_executable(
  InformatiqueGraphique
  main.cpp
  ${SRC_DIR}/Camera/Camera.cpp
  ${SRC_DIR}/Vector/Vector.cpp
  ${SRC_DIR}/Ray/Ray.cpp
  ${SRC_DIR}/Sphere/Sphere.cpp
  ${SRC_DIR}/Environment/Environment.cpp
  ${SRC_DIR}/Mesh/TriangleIndices.cpp
  ${SRC_DIR}/Mesh/TriangleMesh.cpp
  ${SRC_DIR}/Object/Object.cpp
)