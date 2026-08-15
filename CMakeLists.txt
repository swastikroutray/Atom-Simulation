cmake_minimum_required(VERSION 3.15)
project(MyProject)
set(CMAKE_CXX_STANDARD 17)

find_package(glfw3 CONFIG REQUIRED)
find_package(GLEW REQUIRED)
find_package(glm CONFIG REQUIRED)
find_package(OpenGL REQUIRED)

add_executable(MyProject src/atom.cpp)

target_link_libraries(MyProject PRIVATE
    glfw GLEW::GLEW glm::glm OpenGL::GL
)