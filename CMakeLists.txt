cmake_minimum_required(VERSION 3.11) # FetchContent is available in 3.11+
project(DivFree_SPH)

# Generate compile_commands.json
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Dependencies
find_package(raylib QUIET)
if (NOT raylib_FOUND)
    include(FetchContent)
    FetchContent_Declare(
        raylib
        GIT_REPOSITORY https://github.com/raysan5/raylib.git
        GIT_TAG 4.2.0
    )
    FetchContent_MakeAvailable(raylib)
endif()

# raylib-cpp
find_package(raylib_cpp QUIET)
if (NOT raylib_cpp_FOUND)
    include(FetchContent)
    FetchContent_Declare(
        raylib_cpp
        GIT_REPOSITORY https://github.com/RobLoach/raylib-cpp.git
        GIT_TAG v4.2.6
    )
    FetchContent_MakeAvailable(raylib_cpp)
endif()

# Our Project
FILE(GLOB_RECURSE SOURCES "src/*.cpp" "src/*.hpp")
FILE(GLOB_RECURSE HEADERS "src/*.h")
FILE(GLOB_RECURSE GLSL "resources/*.glsl")


add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS} ${GLSL})
#set(raylib_VERBOSE 1)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 20)
target_link_libraries(${PROJECT_NAME} PUBLIC raylib raylib_cpp)

# Web Configurations
if (${PLATFORM} STREQUAL "Web")
    # Tell Emscripten to build an example.html file.
    set_target_properties(${PROJECT_NAME} PROPERTIES SUFFIX ".html")
endif()

# Checks if OSX and links appropriate frameworks (Only required on MacOS)
if (APPLE)
    target_link_libraries(${PROJECT_NAME} "-framework IOKit")
    target_link_libraries(${PROJECT_NAME} "-framework Cocoa")
    target_link_libraries(${PROJECT_NAME} "-framework OpenGL")
endif()