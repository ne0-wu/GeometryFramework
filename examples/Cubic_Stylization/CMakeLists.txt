add_executable(cubic_stylization main.cpp)

target_link_libraries(cubic_stylization GeometryFw)

set_target_properties(cubic_stylization PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/cubic_stylization)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/cubic_stylization)

add_custom_command(
	TARGET cubic_stylization POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E copy_directory
		"${CMAKE_SOURCE_DIR}/shaders"
		"$<TARGET_FILE_DIR:cubic_stylization>/shaders"
	COMMAND ${CMAKE_COMMAND} -E copy_directory
		"${CMAKE_SOURCE_DIR}/meshes"
		"$<TARGET_FILE_DIR:cubic_stylization>/meshes"
)