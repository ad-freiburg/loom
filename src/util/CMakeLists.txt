file(GLOB_RECURSE util_SRC *.cpp)
list(REMOVE_ITEM util_SRC TestMain.cpp)
add_library(util ${util_SRC})

if (!BZIP2_FOUND)
    find_package(BZip2)
    if (ZLIB_FOUND)
		add_definitions( -DBZLIB_FOUND=${BZIP2_FOUND} )
	endif()
endif()

if (BZIP2_FOUND)
	include_directories( ${BZIP2_INCLUDE_DIR} )
	target_link_libraries( util ${BZIP2_LIBRARIES} )
endif(BZIP2_FOUND)

if (!ZLIB_FOUND)
    find_package(ZLIB)
    if (ZLIB_FOUND)
		add_definitions( -DZLIB_FOUND=${ZLIB_FOUND} )
	endif()
endif()

if (ZLIB_FOUND)
	include_directories( ${ZLIB_INCLUDE_DIRS} )
	target_link_libraries( util ${ZLIB_LIBRARIES} )
endif(ZLIB_FOUND)

add_subdirectory(tests)
