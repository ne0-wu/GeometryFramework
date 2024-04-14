cmake_minimum_required(VERSION 3.25)

set(CMAKE_CXX_STANDARD 17)

project(HelloOpenGL)

# Setup libraries
# ---------------
find_package(OpenGL)	# Find the OpenGL package
find_package(OpenMP)	# Find the OpenMP package

include(FetchContent)

# GLFW
FetchContent_Declare(
	glfw
	GIT_REPOSITORY	https://github.com/glfw/glfw.git
	GIT_TAG			3.4
)
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_INSTALL OFF CACHE BOOL "" FORCE)

# GLAD
FetchContent_Declare(
	glad
	GIT_REPOSITORY	https://github.com/Dav1dde/glad
	GIT_TAG			glad2
	SOURCE_SUBDIR	cmake
)

# OpenMesh
FetchContent_Declare(
	openmesh
	GIT_REPOSITORY	https://www.graphics.rwth-aachen.de:9000/OpenMesh/OpenMesh.git
	GIT_TAG			OpenMesh-10.0
)

FetchContent_MakeAvailable(
	glfw
	glad
	openmesh
)

glad_add_library(glad_gl_core_33 REPRODUCIBLE API gl:core=3.3)

add_compile_definitions(_USE_MATH_DEFINES)
add_compile_definitions(OM_STATIC_BUILD)

# Include directories
include_directories(include)
include_directories(3rdparty/include)

# Find source files
file(GLOB_RECURSE SOURCE_FILES src/*.cpp)

add_library(GeometryFw STATIC ${SOURCE_FILES})

# Link the libraries
target_link_libraries(GeometryFw PUBLIC
	${OPENGL_LIBRARIES}	# OpenGL libraries
	glad_gl_core_33		# GLAD library for loading OpenGL functions
	glfw				# GLFW library for creating windows and receiving input and events
	OpenMeshCore		# Core OpenMesh library
	OpenMeshTools		# Tools provided by OpenMesh
	OpenMP::OpenMP_CXX	# OpenMP library for parallel computing
)

add_subdirectory(examples)

# # Handle excutable output path and shaders
# set_property(TARGET HelloOpenGL PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "$<TARGET_FILE_DIR:HelloOpenGL>")

# file(REMOVE_RECURSE ${CMAKE_BINARY_DIR}/HelloOpenGL/)

# set_target_properties(HelloOpenGL PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/HelloOpenGL)

# set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/HelloOpenGL)

# add_custom_command(
# 	TARGET HelloOpenGL POST_BUILD
# 	COMMAND ${CMAKE_COMMAND} -E copy_directory
# 		"${CMAKE_SOURCE_DIR}/shaders"
# 		"$<TARGET_FILE_DIR:HelloOpenGL>/shaders"
# 	COMMAND ${CMAKE_COMMAND} -E copy_directory
# 		"${CMAKE_SOURCE_DIR}/meshes"
# 		"$<TARGET_FILE_DIR:HelloOpenGL>/meshes"
# )