include_directories(
	${TRANSITMAP_INCLUDE_DIR}
	)

add_executable(loomTest TestMain.cpp)
target_link_libraries(loomTest loom_dep shared_dep dot_dep util ${GLPK_LIBRARY} ${GUROBI_LIBRARY} ${COIN_LIBRARIES} -lpthread)
