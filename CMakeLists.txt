cmake_minimum_required (VERSION 2.8)

project (transitmap)

enable_testing()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/build")

find_package(Boost COMPONENTS system program_options thread filesystem REQUIRED)
include_directories("build" ${Boost_INCLUDE_DIRS})

find_package(OpenMP)
if(OPENMP_FOUND)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

find_package(GDAL REQUIRED)
find_package(GLPK REQUIRED)

# set compiler flags, see http://stackoverflow.com/questions/7724569/debug-vs-release-in-cmake
set(CMAKE_CXX_FLAGS                "-fopenmp -DLOGLEVEL=3 -O3 -Wall -Wno-format-extra-args -Wextra -Wformat-nonliteral -Wformat-security -Wformat=2")
set(CMAKE_CXX_FLAGS_DEBUG          "-Og -g -DLOGLEVEL=5")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "-DLOGLEVEL=3")
set(CMAKE_CXX_FLAGS_RELEASE        "-DLOGLEVEL=3")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -DLOGLEVEL=5")

# export compile commands to tools like clang
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Compiler-specific C++11 activation.
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    execute_process(
        COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
    if ((GCC_VERSION VERSION_GREATER 4.7 OR GCC_VERSION VERSION_EQUAL 4.7))
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    elseif (GCC_VERSION VERSION_EQUAL 4.6)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
    else ()
        message(FATAL_ERROR "${PROJECT_NAME} requires g++ 4.6 or greater.")
    endif ()
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++ -std=c++11")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else ()
    message(FATAL_ERROR "Your C++ compiler does not support C++11.")
endif ()

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

add_test(ctest_build_skeletonbuilder_test "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target skeletonbuilderTest)
add_test(skeletonbuilder_test ${EXECUTABLE_OUTPUT_PATH}/skeletonbuilderTest)
set_tests_properties (skeletonbuilder_test PROPERTIES DEPENDS ctest_build_skeletonbuilder_test)

add_test(ctest_build_transitmapper_test "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target transitmapperTest)
add_test(transitmapper_test ${EXECUTABLE_OUTPUT_PATH}/transitmapperTest)
set_tests_properties (transitmapper_test PROPERTIES DEPENDS ctest_build_transitmapper_test)


add_test(ctest_build_pbutil_test "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target pbutilTest)
add_test(pbutil_test ${EXECUTABLE_OUTPUT_PATH}/pbutilTest)
set_tests_properties (pbutil_test PROPERTIES DEPENDS ctest_build_pbutil_test)

# handles install target

install(
  FILES README.md DESTINATION share/${PROJECT_NAME} PERMISSIONS WORLD_READ
)

install(
  FILES build/transitmap DESTINATION bin
  PERMISSIONS WORLD_EXECUTE
)

