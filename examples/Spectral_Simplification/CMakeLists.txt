add_executable(spectral_simp main.cpp)

target_link_libraries(spectral_simp GeometryFw)

set_property(TARGET spectral_simp PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "$<TARGET_FILE_DIR:spectral_simp>")
set_target_properties(spectral_simp PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/spectral_simp)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/spectral_simp)

add_custom_command(
	TARGET spectral_simp POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E copy_directory
		"${CMAKE_SOURCE_DIR}/shaders"
		"$<TARGET_FILE_DIR:spectral_simp>/shaders"
	COMMAND ${CMAKE_COMMAND} -E copy_directory
		"${CMAKE_SOURCE_DIR}/meshes"
		"$<TARGET_FILE_DIR:spectral_simp>/meshes"
)