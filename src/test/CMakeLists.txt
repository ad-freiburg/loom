file(GLOB_RECURSE transitmap_TEST *.cpp)

include_directories(
	${TRANSITMAP_INCLUDE_DIR}
)

add_executable(transitmapTest ${transitmap_TEST})
target_link_libraries(transitmapTest transitmap_dep ${Boost_LIBRARIES} -lproj)
