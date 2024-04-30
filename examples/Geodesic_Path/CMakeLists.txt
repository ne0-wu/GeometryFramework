add_executable(geodesic_path main.cpp)

target_link_libraries(geodesic_path GeometryFw)

set_target_properties(geodesic_path PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/geodesic_path)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/geodesic_path)

add_custom_command(
	TARGET geodesic_path POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E copy_directory
		"${CMAKE_SOURCE_DIR}/shaders"
		"$<TARGET_FILE_DIR:geodesic_path>/shaders"
	COMMAND ${CMAKE_COMMAND} -E copy_directory
		"${CMAKE_SOURCE_DIR}/meshes"
		"$<TARGET_FILE_DIR:geodesic_path>/meshes"
)