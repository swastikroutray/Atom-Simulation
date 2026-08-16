# Atom Simulation

***It is a real-time hydrogen quantum orbital visualizer written in C++ and OpenGL. Electron probability density clouds are generated directly from the hydrogen wavefunction with no pre-baked meshes, just physics sampled live.***

<img width="1914" height="1075" alt="image" src="https://github.com/user-attachments/assets/100ac5de-2d9a-42ea-a517-ace30a9f45ce" />


## How it works

***The visualizer solves the hydrogen atom's wavefunction and samples particle positions from the resulting probability density:***

**i) Radial distribution** — ***computed from the Associated Laguerre polynomials, sampled via a cumulative distribution function (CDF)***

**ii) Angular distribution** — ***computed from the Associated Legendre polynomials (spherical harmonics), also sampled via CDF***

**iii) Particle** - ***Each particle is colored by local probability density using a heatmap gradient***

**iv) Plane** - ***Three reference planes (XY, XZ, YZ) are rendered around the cloud, scaling automatically with the size of the orbital***

***Every time the quantum numbers change, the probability distributions are recalculated, and a fresh set of particles is sampled; so the shape you see is always mathematically accurate for the current n, l, m.***

## Controls

***Key	Action***

**W / S**	               ***Increase/decrease n (principal quantum number)***

**E / D**	               ***Increase/decrease l (azimuthal quantum number)***

**R / F**	               ***Increase/decrease m (magnetic quantum number)***

**Left-click + drag**	   ***Orbit the camera around the atom***

**Scroll**	             ***Zoom in/out***

**Quantum numbers are automatically clamped to valid combinations (0 ≤ l < n, -l ≤ m ≤ l).**

## Tech stack

***C++17***
***OpenGL (core profile) for rendering***
***GLFW for windowing and input***
***GLEW for OpenGL function loading***
***GLM for vector/matrix math***
***CMake + vcpkg for building and dependency management***

## Building from source

***This repo tracks only the source files in src/. To build it locally, create the following two files in the project root:***

***vcpkg.json***

***json
{
  "name": "atom-simulation",
  "version": "1.0.0",
  "dependencies": ["glfw3", "glew", "glm"]
}***

***CMakeLists.txt***

***cmake
cmake_minimum_required(VERSION 3.15)
project(Atoms CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
find_package(OpenGL REQUIRED)
find_package(GLEW REQUIRED)
find_package(glfw3 REQUIRED)
find_package(glm REQUIRED)
include_directories(src ${OPENGL_INCLUDE_DIRS} ${GLEW_INCLUDE_DIRS})
add_executable(atom_realtime src/atom_realtime.cpp)
target_link_libraries(atom_realtime PRIVATE OpenGL::GL GLEW::GLEW glfw glm::glm)
add_executable(atom_raytracer src/atom_raytracer.cpp)
target_link_libraries(atom_raytracer PRIVATE OpenGL::GL GLEW::GLEW glfw glm::glm)***

***Then, with vcpkg installed:***

**bash**

**vcpkg install**

**cmake -B build -S . -DCMAKE_TOOLCHAIN_FILE=<path-to-vcpkg>/scripts/buildsystems/vcpkg.cmake**

**cmake --build build**

**The compiled executable will be in build/.**

## Source files

**File	Description**

***src/atom_realtime.cpp	Real-time 3D orbital visualizer with orbit camera and reference planes***

***src/atom_raytracer.cpp	Raytraced rendering of the orbital***

***src/atom.cpp	2D Bohr model visualization***

***src/wave_atom.cpp	2D electron wave visualization***
