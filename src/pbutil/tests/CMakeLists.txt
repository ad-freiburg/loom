include_directories(
	${TRANSITMAP_INCLUDE_DIR}
)

add_executable(pbutilTest TestMain.cpp)
target_link_libraries(pbutilTest)
