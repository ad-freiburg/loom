include_directories(
	${LOOM_INCLUDE_DIR}
)

add_executable(octiTest TestMain.cpp)
target_link_libraries(octiTest octi_dep util ad_cppgtfs)
