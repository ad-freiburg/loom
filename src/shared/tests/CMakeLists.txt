include_directories(
	${LOOM_INCLUDE_DIR}
)

add_executable(sharedTest TestMain.cpp)

target_link_libraries(sharedTest shared_dep util ${GUROBI_LIBRARY} ${GLPK_LIBRARY} ${COIN_LIBRARIES})
