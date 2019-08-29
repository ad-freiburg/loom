file(GLOB_RECURSE test_SRC *.cpp)
list(REMOVE_ITEM test_SRC TestMain.cpp)

include_directories(
	${TRANSITMAP_INCLUDE_DIR}
)

add_executable(topoTest TestMain.cpp)
add_library(topo_test_dep ${test_SRC})
target_link_libraries(topoTest topo_test_dep topo_dep shared_dep dot_dep transitmap_dep util ${Boost_LIBRARIES} -lproj)
