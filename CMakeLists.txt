# SPDX-License-Identifier: BSD-2-Clause

# Copyright (c) 2021 NKI/AVL, Netherlands Cancer Institute

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

cmake_minimum_required(VERSION 3.15)

# set the project name
project(alphafill VERSION 1.0.0 LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(GNUInstallDirs)
include(CheckFunctionExists)
include(CheckIncludeFiles)
include(CheckLibraryExists)
include(CMakePackageConfigHelpers)
include(Dart)
include(FindFilesystem)
include(FindPkgConfig)
include(GenerateExportHeader)

set(CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Filesystem REQUIRED)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-unused-parameter -Wno-missing-field-initializers")
elseif(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
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

include(GetGitRevisionDescription)
option(GENERATE_CUSTOM_VERSION "Generate a custom version string")

if(NOT GENERATE_CUSTOM_VERSION OR GIT-NOTFOUND OR HEAD-HASH-NOTFOUND)
	SET(BUILD_VERSION_STRING ${PROJECT_VERSION})
else()
	get_git_head_revision(REFSPEC COMMITHASH)

	# Generate our own version string
	git_describe_working_tree(BUILD_VERSION_STRING --match=build --dirty)
endif()

# generate version.h
string(TIMESTAMP BUILD_DATE_TIME "%Y-%m-%dT%H:%M:%SZ" UTC)
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
	add_compile_definitions(USE_RSRC WEBAPP_USES_RESOURCES)

	mrc_write_header(${CMAKE_BINARY_DIR}/mrsrc.hpp)
endif()

set(CMAKE_THREAD_PREFER_PTHREAD)
set(THREADS_PREFER_PTHREAD_FLAG)
find_package(Threads REQUIRED)

set(Boost_DETAILED_FAILURE_MSG ON)

find_package(zeep 5.1.6 REQUIRED)
find_package(cifpp REQUIRED)
find_package(Boost COMPONENTS system date_time program_options regex iostreams)

# The alphafill executable

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

target_include_directories(alphafill PRIVATE ${CMAKE_BINARY_DIR})

target_link_libraries(alphafill cifpp::cifpp ZLIB::ZLIB BZip2::BZip2 zeep::zeep ${CMAKE_THREAD_LIBS_INIT})

if(USE_RSRC)
	mrc_target_resources(alphafill ${CIFPP_SHARE_DIR}/mmcif_pdbx_v50.dic)
endif()

# The web interface, af-filledd

pkg_check_modules(PQ libpq REQUIRED)
pkg_check_modules(PQXX libpqxx>=7.2 REQUIRED)
# pkg_check_modules(MAILIO mailio REQUIRED)
# pkg_check_modules(SSL libssl)
# pkg_check_modules(CRYPTO libcrypto)

find_program(YARN yarn REQUIRED)

pkg_check_modules(ZLIB zlib)
# pkg_check_modules(BZIP2 libbz2)
find_library(BZIP2_LIBRARIES bz2)

find_program(YARN yarn REQUIRED)

# yarn rules for javascripts

set(webpack_input
	${CMAKE_SOURCE_DIR}/webapp/index.js
	${CMAKE_SOURCE_DIR}/webapp/model.js
	${CMAKE_SOURCE_DIR}/webapp/lists.js
	${CMAKE_SOURCE_DIR}/webapp/molstar.tsx
	${CMAKE_SOURCE_DIR}/scss/pdb-redo-bootstrap.scss
)

set(webpack_output ${CMAKE_SOURCE_DIR}/docroot/scripts/index.js)

add_custom_command(
	OUTPUT ${webpack_output}
	BYPRODUCTS 
	DEPENDS ${webpack_input}
	COMMAND ${YARN} run $<IF:$<CONFIG:Debug>,build,build-production>
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

add_custom_target(webpack DEPENDS ${webpack_output})

add_executable(af-filledd
	${CMAKE_SOURCE_DIR}/src/af-filledd.cpp
	${CMAKE_SOURCE_DIR}/src/db-connection.cpp
	${CMAKE_SOURCE_DIR}/src/data-service.cpp
)

target_include_directories(af-filledd PRIVATE zeep::zeep cifpp::cifpp squeeze::squeeze ${CMAKE_SOURCE_DIR}/include ${Boost_INCLUDE_DIRS} ${CMAKE_BINARY_DIR})
target_link_libraries(af-filledd
	PRIVATE zeep::zeep cifpp::cifpp Boost::date_time Boost::program_options Boost::regex Boost::iostreams
	${PQXX_LIBRARIES} ${PQ_LIBRARIES} ${MAILIO_LINK_LIBRARIES} ${SSL_LIBRARIES} ${CRYPTO_LIBRARIES} ${ZLIB_LIBRARIES} ${BZIP2_LIBRARIES} std::filesystem
	${CMAKE_THREAD_LIBS_INIT})

add_dependencies(af-filledd webpack)

if(USE_RSRC)
	mrc_target_resources(af-filledd ${CMAKE_SOURCE_DIR}/docroot/ ${CMAKE_SOURCE_DIR}/db-schema.sql)
endif()

install(TARGETS alphafill af-filledd
    RUNTIME DESTINATION ${BIN_INSTALL_DIR}
)
