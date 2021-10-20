cmake_minimum_required(VERSION 3.15)

# set the project name
project(alphafill VERSION 1.0.0 LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(GNUInstallDirs)
include(Dart)
include(FindFilesystem)

set(CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Filesystem REQUIRED)

if(CMAKE_COMPILER_IS_GNUCC)
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-unused-parameter")
endif()
if(MSVC)
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} /W4")
endif()

if(NOT "$ENV{CCP4}" STREQUAL "")
	set(CCP4 $ENV{CCP4})
	list(PREPEND CMAKE_MODULE_PATH "${CCP4}/Lib")
	list(APPEND CMAKE_PREFIX_PATH ${CCP4})
	if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
		set(CMAKE_PREFIX_PATH ${CCP4})
	endif()
endif()

if(MSVC)
    # make msvc standards compliant...
    add_compile_options(/permissive-)

	macro(get_WIN32_WINNT version)
		if (WIN32 AND CMAKE_SYSTEM_VERSION)
			set(ver ${CMAKE_SYSTEM_VERSION})
			string(REPLACE "." "" ver ${ver})
			string(REGEX REPLACE "([0-9])" "0\\1" ver ${ver})

			set(${version} "0x${ver}")
		endif()
	endmacro()

	get_WIN32_WINNT(ver)
	add_definitions(-D_WIN32_WINNT=${ver})

	# On Windows, do not install in the system location
	if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT AND NOT ($ENV{LOCALAPPDATA} STREQUAL ""))
		message(STATUS "The executable will be installed in $ENV{LOCALAPPDATA}")
		set(CMAKE_INSTALL_PREFIX "$ENV{LOCALAPPDATA}" CACHE PATH "..." FORCE)
	endif()
endif()

if(UNIX AND NOT APPLE)
	# On Linux, install in the $HOME/.local folder by default
	if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
		message(STATUS "The library and auxiliary files will be installed in $ENV{HOME}/.local")
		set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/.local" CACHE PATH "..." FORCE)
	endif()
endif()

# Create a revision file, containing the current git version info

find_package(Git)
if(GIT_FOUND AND EXISTS "${CMAKE_SOURCE_DIR}/.git")
	include(GetGitRevisionDescription)
    get_git_head_revision(REFSPEC COMMITHASH)

	# Generate our own version string
	git_describe_working_tree(BUILD_VERSION_STRING --match=build --dirty)
else()
    message(STATUS "Git not found, cannot set version info")

    SET(BUILD_VERSION_STRING ${PROJECT_VERSION})
endif()

# generate version.h
include_directories(${CMAKE_BINARY_DIR} PRIVATE)
string(TIMESTAMP BUILD_DATE_TIME "%Y-%m-%d" UTC)
configure_file("${CMAKE_SOURCE_DIR}/src/revision.hpp.in" "${CMAKE_BINARY_DIR}/revision.hpp" @ONLY)

# Optionally use mrc to create resources
find_package(Mrc)

if(MRC_FOUND)
	option(USE_RSRC "Use mrc to create resources" ON)
else()
	message(STATUS "Not using resources since mrc was not found")
endif()

if(USE_RSRC)
	set(USE_RSRC 1)

	message("Using resources compiled with ${MRC_EXECUTABLE}")
	add_compile_definitions(USE_RSRC)
endif()

set(CMAKE_THREAD_PREFER_PTHREAD)
set(THREADS_PREFER_PTHREAD_FLAG)
find_package(Threads REQUIRED)

set (Boost_DETAILED_FAILURE_MSG ON)
find_package(Boost 1.70.0 REQUIRED COMPONENTS program_options system iostreams date_time regex)

find_package(ZLIB)
find_package(BZip2)
find_package(cifpp 3.0.0 CONFIG REQUIRED HINTS $ENV{LOCALAPPDATA}/cifpp)
find_package(zeep 5.1.4 REQUIRED)

if(CIFPP_FOUND)
	add_compile_definitions("DATA_DIR=\"${CIFPP_SHARE_DIR}\"")
endif()

add_executable(alphafill
	${PROJECT_SOURCE_DIR}/src/blast.cpp
	${PROJECT_SOURCE_DIR}/src/matrix.cpp
	${PROJECT_SOURCE_DIR}/src/alphafill.cpp)

target_include_directories(alphafill
	PUBLIC
	cifpp::cifpp
	zeep::zeep
	${CMAKE_SOURCE_DIR}/include
	${cifpp_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})

target_link_libraries(alphafill cifpp::cifpp ZLIB::ZLIB BZip2::BZip2 zeep::zeep ${CMAKE_THREAD_LIBS_INIT})

if(USE_RSRC)
	mrc_target_resources(alphafill ${CIFPP_SHARE_DIR}/mmcif_pdbx_v50.dic)
endif()

install(TARGETS ${PROJECT_NAME}
    RUNTIME DESTINATION ${BIN_INSTALL_DIR}
)

# # manual

# if(UNIX)
# 	install(FILES doc/alphafill.1
# 		DESTINATION ${CMAKE_INSTALL_DATADIR}/man/man1)
# endif()

# if(EXISTS "${CCP4}/html")
# 	install(FILES doc/alphafill.html
# 		DESTINATION ${CCP4}/html)
# endif()

# test

add_executable(unit-test ${PROJECT_SOURCE_DIR}/test/unit-test.cpp)

target_include_directories(unit-test PRIVATE
	${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(unit-test Threads::Threads ${Boost_LIBRARIES} cifpp::cifpp)

if(${ZLIB_FOUND})
	target_link_libraries(unit-test ZLIB::ZLIB)
endif()

if(${BZip2_FOUND})
	target_link_libraries(unit-test BZip2::BZip2)
endif()

if(MSVC)
	# Specify unwind semantics so that MSVC knowns how to handle exceptions
	target_compile_options(unit-test PRIVATE /EHsc)
endif()

enable_testing()

add_test(NAME unit-test
	COMMAND $<TARGET_FILE:unit-test>
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test)