cmake_minimum_required(VERSION 2.8)
project(runThis)

set(SOURCES
    src/fluidSimMain3.cpp
    src/Classes/Shader.cpp
    src/Classes/Texture.cpp
    src/Particles/ParticleSystem.cpp
    src/Particles/Particle.cpp
    src/MarchingCubes/MarchingCubes.cpp
    src/Classes/TransformFeedback.cpp
    src/Classes/Rainbow.cpp
    src/FluidSim/Image.cpp
)
set(HEADERS
    src/Classes/Shader.hpp
    src/Classes/Texture.hpp
    src/Classes/Camera.hpp
    src/Particles/ParticleSystem.hpp
    src/Particles/Particle.hpp
    src/MarchingCubes/MarchingCubes.hpp
    src/constants.hpp
    src/Classes/TransformFeedback.hpp
    src/Classes/Rainbow.hpp
    src/FluidSim/Image.h
    src/FluidSim/types.h
)
set(SHADERS
)
set(RESOURCES
)
set(CMAKE_CXX_STANDARD 11)

add_executable(${CMAKE_PROJECT_NAME} ${SOURCES} ${HEADERS})
add_custom_target(run
    COMMAND ${CMAKE_BINARY_DIR}/${CMAKE_PROJECT_NAME}
    DEPENDS ${CMAKE_BINARY_DIR}/${CMAKE_PROJECT_NAME}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

include_directories(${CMAKE_PROJECT_NAME} ${CMAKE_SOURCE_DIR}/src)
include_directories(${CMAKE_PROJECT_NAME} ${CMAKE_SOURCE_DIR}/src/imgui)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "/Users/ashleydattalo/cmake/modules/")
find_package(GLFW REQUIRED)
if (GLFW_FOUND)
    include_directories(${CMAKE_PROJECT_NAME} ${GLFW_INCLUDE_DIR})
    target_link_libraries(${CMAKE_PROJECT_NAME} ${GLFW_LIBRARY})
else()
    message(FATAL_ERROR "Can't find GLFW")
endif()

find_package(OpenGL REQUIRED)
if (OPENGL_FOUND)
    include_directories(${CMAKE_PROJECT_NAME} ${OPENGL_INCLUDE_DIR})
    target_link_libraries(${CMAKE_PROJECT_NAME} ${OPENGL_gl_LIBRARY})
else()
    message(FATAL_ERROR "Can't find OpenGL")
endif()

find_package(GLEW)
if (GLEW_FOUND)
    include_directories(${CMAKE_PROJECT_NAME} ${GLEW_INCLUDE_DIR})
    target_link_libraries(${CMAKE_PROJECT_NAME} ${GLEW_LIBRARIES})
elseif (DEFINED ENV{GLEW_INCLUDE_DIR})
    message(STATUS "using GLEW_INCLUDE_DIR environment variable")
    include_directories(${CMAKE_PROJECT_NAME} "$ENV{GLEW_INCLUDE_DIR}/include")
    
    if (APPLE)
       target_link_libraries(${CMAKE_PROJECT_NAME} "$ENV{GLEW_INCLUDE_DIR}/lib/libGLEW.a")
    elseif (WIN32)
       target_link_libraries(${CMAKE_PROJECT_NAME} "$ENV{GLEW_INCLUDE_DIR}/lib/Release/Win32/glew32s.lib")
    endif()
else()
    message(FATAL_ERROR "Can't find GLEW")
endif()

if (DEFINED ENV{GLM_LOCATION})
    message(STATUS "using GLM_LOCATION environment variable")
    include_directories(${CMAKE_PROJECT_NAME} "$ENV{GLM_LOCATION}")
endif()

find_package(SFML COMPONENTS audio system)
if (SFML_FOUND)
    include_directories(${CMAKE_PROJECT_NAME} ${SFML_INCLUDE_DIR})
    target_link_libraries(${CMAKE_PROJECT_NAME} ${SFML_LIBRARIES})
else()
    message(FATAL_ERROR "Can't find SFML. Set environment variable SFML_ROOT.")
endif()
