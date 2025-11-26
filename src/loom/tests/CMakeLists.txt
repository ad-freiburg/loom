include_directories(
	${LOOM_INCLUDE_DIR}
)

add_executable(loomTest TestMain.cpp)
target_link_libraries(loomTest loom_dep shared_dep dot_dep pb_util_graph pb_util pb_util_json pb_util_geo ${GLPK_LIBRARY} ${GUROBI_LIBRARY} ${COIN_LIBRARIES} -lpthread)
