cmake_minimum_required(VERSION 3.25)

set(CMAKE_CXX_STANDARD 17)

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")

project(HelloOpenGL)

find_package(OpenGL) # Find the OpenGL package
find_package(OpenMP) # Find the OpenMP package

# Setup libraries
include(FetchContent)

# GLFW
FetchContent_Declare(
	glfw
	GIT_REPOSITORY	https://github.com/glfw/glfw.git
	GIT_TAG			3.4
)
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_INSTALL OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(glfw)

# GLAD
FetchContent_Declare(
	glad
	GIT_REPOSITORY	https://github.com/Dav1dde/glad
	GIT_TAG			glad2
	SOURCE_SUBDIR	cmake
)
FetchContent_MakeAvailable(glad)
glad_add_library(glad_gl_core_33 REPRODUCIBLE API gl:core=3.3)

# OpenMesh
FetchContent_Declare(
	openmesh
	GIT_REPOSITORY	https://www.graphics.rwth-aachen.de:9000/OpenMesh/OpenMesh.git
	GIT_TAG			OpenMesh-10.0
)
if (WIN32 AND "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wa,-mbig-obj")
endif()
FetchContent_MakeAvailable(openmesh)

add_compile_definitions(_USE_MATH_DEFINES)
add_compile_definitions(OM_STATIC_BUILD)

# Include directories
include_directories(include)
include_directories(3rdparty/include)

# Find source files
file(GLOB_RECURSE SOURCE_FILES src/*.cpp)

# Add the executable
add_executable(HelloOpenGL
	${SOURCE_FILES}
	# 3rdparty/src/glad.c
)

# Link the libraries
target_link_libraries(HelloOpenGL
	${OPENGL_LIBRARIES}	# OpenGL libraries
	glad_gl_core_33		# GLAD library for loading OpenGL functions
	OpenMP::OpenMP_CXX	# OpenMP
	glfw				# GLFW library for creating windows, contexts and surfaces, receiving input and events
	OpenMeshCore		# Core OpenMesh library
	OpenMeshTools		# Tools provided by OpenMesh
)

# Handle excutable output path and shaders
set_property(TARGET HelloOpenGL PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "$<TARGET_FILE_DIR:HelloOpenGL>")

file(REMOVE_RECURSE ${CMAKE_BINARY_DIR}/HelloOpenGL/)

set_target_properties(HelloOpenGL PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/HelloOpenGL)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/HelloOpenGL)

add_custom_command(
	TARGET HelloOpenGL POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E copy_directory
		"${CMAKE_SOURCE_DIR}/shaders"
		"$<TARGET_FILE_DIR:HelloOpenGL>/shaders"
	COMMAND ${CMAKE_COMMAND} -E copy_directory
		"${CMAKE_SOURCE_DIR}/meshes"
		"$<TARGET_FILE_DIR:HelloOpenGL>/meshes"
)