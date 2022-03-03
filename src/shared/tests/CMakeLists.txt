include_directories(
	${TRANSITMAP_INCLUDE_DIR}
)

add_executable(sharedTest TestMain.cpp)

target_link_libraries(sharedTest shared_dep util ${Boost_LIBRARIES} ${GUROBI_LIBRARY} ${GLPK_LIBRARY} ${COIN_LIBRARIES})
