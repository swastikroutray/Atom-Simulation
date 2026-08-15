cmake_minimum_required(VERSION 3.15)
project(Atoms)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(glfw3 CONFIG REQUIRED)
find_package(GLEW REQUIRED)
find_package(glm CONFIG REQUIRED)
find_package(OpenGL REQUIRED)

add_executable(Atoms src/atom.cpp)

target_link_libraries(Atoms PRIVATE
    glfw
    GLEW::GLEW
    glm::glm
    OpenGL::GL
)