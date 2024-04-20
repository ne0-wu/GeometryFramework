cmake_minimum_required(VERSION 3.25)

set(CMAKE_CXX_STANDARD 17)

project(GeometryFramework)

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

# Dear ImGui
FetchContent_Declare(
	imgui
	GIT_REPOSITORY	https://github.com/ocornut/imgui
	GIT_TAG			v1.90.5-docking
)

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
set(BUILD_APPS OFF CACHE BOOL "" FORCE)
# set(OPENMESH_DOCS OFF CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(
	glfw
	imgui
	glad
	openmesh
)

add_library(imgui STATIC
	${imgui_SOURCE_DIR}/imgui.cpp
	${imgui_SOURCE_DIR}/imgui_draw.cpp
	${imgui_SOURCE_DIR}/imgui_widgets.cpp
	${imgui_SOURCE_DIR}/imgui_demo.cpp
	${imgui_SOURCE_DIR}/imgui_tables.cpp
	${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp
	${imgui_SOURCE_DIR}/backends/imgui_impl_opengl3.cpp
)
target_link_libraries(imgui PUBLIC glfw glad_gl_core_33)
target_include_directories(imgui PUBLIC ${imgui_SOURCE_DIR} ${imgui_SOURCE_DIR}/backends)

glad_add_library(glad_gl_core_33 REPRODUCIBLE API gl:core=3.3)

# Compile definitions for OpenMesh
add_compile_definitions(_USE_MATH_DEFINES)
add_compile_definitions(OM_STATIC_BUILD)

# ==================================================
# Compile definitions for Geometry Processing
# Example: define IMPLEMENT_QEM to use Quadric Error Metrics
# --------------------------------------------------
# add_compile_definitions(IMPLEMENT_QEM)
# add_compile_definitions(IMPLEMENT_PSR)
# add_compile_definitions(IMPLEMENT_LGP)
# ==================================================

# Include directories
include_directories(include)
include_directories(3rdparty/include)

# Compile the source files
file(GLOB_RECURSE SOURCE_FILES src/*.cpp)
add_library(GeometryFw STATIC ${SOURCE_FILES})
target_link_libraries(GeometryFw PUBLIC
	${OPENGL_LIBRARIES}	# OpenGL libraries
	glad_gl_core_33		# GLAD library for loading OpenGL functions
	glfw				# GLFW library for creating windows and receiving input and events
	OpenMeshCore		# Core OpenMesh library
	OpenMeshTools		# Tools provided by OpenMesh
	OpenMP::OpenMP_CXX	# OpenMP library for parallel computing
	imgui				# Dear ImGui library for GUI
)

add_subdirectory(examples)