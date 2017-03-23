include_directories(
	${TRANSITMAP_INCLUDE_DIR}
)

add_executable(transitmapperTest TestMain.cpp)
target_link_libraries(transitmapperTest transitmap_dep ${Boost_LIBRARIES} -lproj)
