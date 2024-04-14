add_executable(local_global main.cpp)

target_link_libraries(local_global GeometryFw)

set_property(TARGET local_global PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "$<TARGET_FILE_DIR:local_global>")
set_target_properties(local_global PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/local_global)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/local_global)

add_custom_command(
	TARGET local_global POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E copy_directory
		"${CMAKE_SOURCE_DIR}/shaders"
		"$<TARGET_FILE_DIR:local_global>/shaders"
	COMMAND ${CMAKE_COMMAND} -E copy_directory
		"${CMAKE_SOURCE_DIR}/meshes"
		"$<TARGET_FILE_DIR:local_global>/meshes"
)