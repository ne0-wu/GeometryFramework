add_executable(qem main.cpp)

target_link_libraries(qem GeometryFw)

set_property(TARGET qem PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "$<TARGET_FILE_DIR:qem>")
set_target_properties(qem PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/qem)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/qem)

add_custom_command(
	TARGET qem POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E copy_directory
		"${CMAKE_SOURCE_DIR}/shaders"
		"$<TARGET_FILE_DIR:qem>/shaders"
	COMMAND ${CMAKE_COMMAND} -E copy_directory
		"${CMAKE_SOURCE_DIR}/meshes"
		"$<TARGET_FILE_DIR:qem>/meshes"
)