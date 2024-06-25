cmake_minimum_required (VERSION 3.5)
set(CMAKE_CXX_STANDARD 11)

project (loom)

if (CMAKE_BUILD_TYPE)
	string(SUBSTRING ${CMAKE_BUILD_TYPE} 0 1 FIRST_CHAR)
	string(TOUPPER ${FIRST_CHAR} FIRST_CHAR)
	string(REGEX REPLACE "^.(.*)" "${FIRST_CHAR}\\1" CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}")
endif()

enable_testing()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}")

find_package(OpenMP)
if(OPENMP_FOUND)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

find_package(GLPK)
find_package(Gurobi)
find_package(COIN)
find_package(Protobuf)
find_package(LibZip)

# set compiler flags, see http://stackoverflow.com/questions/7724569/debug-vs-release-in-cmake
if(OPENMP_FOUND)
	set(CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS} -O3 -Wall -Wno-format-extra-args -Wextra -Wformat-nonliteral -Wno-keyword-macro -Wformat-security -Wformat=2")

else()
	message(WARNING "Configuring without OpenMP!")
	set(CMAKE_CXX_FLAGS "-O3 -Wall -Wno-format-extra-args -Wextra -Wformat-nonliteral -Wno-keyword-macro -Wformat-security -Wformat=2")
endif()

if(GUROBI_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGUROBI_FOUND=1")
endif()

if(Protobuf_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DPROTOBUF_FOUND=1")
endif()

if(GLPK_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGLPK_FOUND=1")
endif()

if(COIN_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCOIN_FOUND=1")
endif()

if (LIBZIP_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DLIBZIP_FOUND=1")
endif()

set(CMAKE_CXX_FLAGS_DEBUG          "-Og -g -DLOGLEVEL=3")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "${CMAKE_CXX_FLAGS} -DLOGLEVEL=2 -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS} -DLOGLEVEL=2 -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS} -g -DLOGLEVEL=3")

# export compile commands to tools like clang
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# http://brianmilco.blogspot.de/2012/11/cmake-automatically-use-git-tags-as.html
include(GetGitRevisionDescription)
git_get_tag(VERSION_GIT)
get_git_is_dirty(VERSION_GIT_IS_DIRTY)
if ("${VERSION_GIT_IS_DIRTY}" STREQUAL "")
    set(VERSION_GIT_FULL "${VERSION_GIT}")
else()
    set(VERSION_GIT_FULL "${VERSION_GIT}-${VERSION_GIT_IS_DIRTY}")
endif()

add_subdirectory(src)

add_test(ctest_build_topo_test "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target topoTest)
add_test(topo_test ${EXECUTABLE_OUTPUT_PATH}/topoTest)
set_tests_properties (topo_test PROPERTIES DEPENDS ctest_build_topo_test)

add_test(ctest_build_util_test "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target utilTest)
add_test(util_test ${EXECUTABLE_OUTPUT_PATH}/utilTest)
set_tests_properties (util_test PROPERTIES DEPENDS ctest_build_util_test)

add_test(ctest_build_shared_test "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target sharedTest)
add_test(shared_test ${EXECUTABLE_OUTPUT_PATH}/sharedTest)
set_tests_properties (shared_test PROPERTIES DEPENDS ctest_build_shared_test)

add_test(ctest_build_loom_test "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target loomTest)
add_test(loom_test ${EXECUTABLE_OUTPUT_PATH}/loomTest)
set_tests_properties (loom_test PROPERTIES DEPENDS ctest_build_loom_test)

add_test(ctest_build_octi_test "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target octiTest)
add_test(octi_test ${EXECUTABLE_OUTPUT_PATH}/octiTest)
set_tests_properties (octi_test PROPERTIES DEPENDS ctest_build_octi_test)

add_test(ctest_build_transitmap_test "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target transitmapTest)
add_test(transitmap_test ${EXECUTABLE_OUTPUT_PATH}/transitmapTest)
set_tests_properties (transitmap_test PROPERTIES DEPENDS ctest_build_transitmap_test)

# handles install target

install(
	FILES README.md DESTINATION share/${PROJECT_NAME} PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
)

install(
  FILES ${CMAKE_BINARY_DIR}/transitmap ${CMAKE_BINARY_DIR}/topo ${CMAKE_BINARY_DIR}/topoeval ${CMAKE_BINARY_DIR}/gtfs2graph ${CMAKE_BINARY_DIR}/loom ${CMAKE_BINARY_DIR}/octi DESTINATION bin
  PERMISSIONS OWNER_EXECUTE GROUP_EXECUTE WORLD_EXECUTE
)

