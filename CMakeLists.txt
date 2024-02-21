# SPDX-License-Identifier: BSD-2-Clause

# Copyright (c) 2021 NKI/AVL, Netherlands Cancer Institute

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# 1. Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimer
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

cmake_minimum_required(VERSION 3.21)

# set the project name
project(
  alphafill
  VERSION 2.1.1
  LANGUAGES C CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(GNUInstallDirs)
include(CheckSymbolExists)
include(CheckFunctionExists)
include(CheckIncludeFiles)
include(CheckLibraryExists)
include(CMakePackageConfigHelpers)
include(GenerateExportHeader)
include(CTest)
include(FetchContent)

set(CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(BUILD_DOCUMENTATION "Build the documentation" OFF)
option(BUILD_WEB_APPLICATION "Build the web application" OFF)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	set(CMAKE_CXX_FLAGS
		"${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-unused-parameter -Wno-missing-field-initializers"
	)
elseif(MSVC)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
endif()

if(WIN32)
	if(${CMAKE_SYSTEM_VERSION} GREATER_EQUAL 10) # Windows 10
		add_definitions(-D _WIN32_WINNT=0x0A00)
	elseif(${CMAKE_SYSTEM_VERSION} EQUAL 6.3) # Windows 8.1
		add_definitions(-D _WIN32_WINNT=0x0603)
	elseif(${CMAKE_SYSTEM_VERSION} EQUAL 6.2) # Windows 8
		add_definitions(-D _WIN32_WINNT=0x0602)
	elseif(${CMAKE_SYSTEM_VERSION} EQUAL 6.1) # Windows 7
		add_definitions(-D _WIN32_WINNT=0x0601)
	elseif(${CMAKE_SYSTEM_VERSION} EQUAL 6.0) # Windows Vista
		add_definitions(-D _WIN32_WINNT=0x0600)
	else() # Windows XP (5.1)
		add_definitions(-D _WIN32_WINNT=0x0501)
	endif()

	add_definitions(-DNOMINMAX)

	# We do not want to write an export file for all our symbols...
	set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

if(MSVC)
	# make msvc standards compliant...
	add_compile_options(/permissive- /bigobj)
	add_link_options(/NODEFAULTLIB:library)

	# static runtime
	set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()

# Create a revision file, containing the current git version info
include(VersionString)
write_version_header("${PROJECT_SOURCE_DIR}/src/")

# Optionally use mrc to create resources
find_package(Mrc QUIET)

if(MRC_FOUND)
	option(USE_RSRC "Use mrc to create resources" ON)
else()
	message(STATUS "Not using resources since mrc was not found")
endif()

# push the BUILD_TESTING flag
set(SAVE_BUILD_TESTING ${BUILD_TESTING})
set(BUILD_TESTING OFF)

find_package(cifpp 7 QUIET)

if(NOT cifpp_FOUND)
  FetchContent_Declare(
    cifpp
    GIT_REPOSITORY https://github.com/PDB-REDO/libcifpp.git
    GIT_TAG v7.0.0)
  FetchContent_MakeAvailable(cifpp)

	set(CIFPP_SHARE_DIR ${cifpp_SOURCE_DIR}/rsrc)
endif()

find_package(libmcfp 1.2.4 QUIET)

if(NOT libmcfp_FOUND)
	FetchContent_Declare(
		libmcfp
		GIT_REPOSITORY https://github.com/mhekkel/libmcfp.git
		GIT_TAG v1.2.4)
	FetchContent_MakeAvailable(libmcfp)
endif()

find_package(zeep 6.0.11 QUIET)

if(NOT zeep_FOUND)
	FetchContent_Declare(
		zeep
		GIT_REPOSITORY https://github.com/mhekkel/libzeep.git
		GIT_TAG v6.0.11)
	FetchContent_MakeAvailable(zeep)
endif()

if(BUILD_WEB_APPLICATION)
	find_package(libpqxx 7.8.0 QUIET)

	if(NOT libpqxx_FOUND)
		FetchContent_Declare(
			libpqxx
			GIT_REPOSITORY https://github.com/jtv/libpqxx
			GIT_TAG 7.8.0)
		FetchContent_MakeAvailable(libpqxx)
	endif()

	find_package(OpenSSL REQUIRED)

	find_program(YARN yarn REQUIRED)

	# Make sure yarn is initialised
	add_custom_command(
		OUTPUT ${PROJECT_SOURCE_DIR}/node_modules
		COMMAND ${YARN}
		WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})

	# yarn rules for javascripts
	set(webpack_input
		${PROJECT_SOURCE_DIR}/webapp/index.js
		${PROJECT_SOURCE_DIR}/webapp/model.js
		${PROJECT_SOURCE_DIR}/webapp/optimized.js
		${PROJECT_SOURCE_DIR}/webapp/lists.js
		${PROJECT_SOURCE_DIR}/webapp/wait.js
		${PROJECT_SOURCE_DIR}/webapp/molstar.tsx
		${PROJECT_SOURCE_DIR}/scss/pdb-redo-bootstrap.scss)

	set(webpack_output ${PROJECT_SOURCE_DIR}/docroot/scripts/index.js)

	add_custom_command(
		OUTPUT ${webpack_output}
		BYPRODUCTS
		DEPENDS ${webpack_input}
		COMMAND ${YARN} run $<IF:$<CONFIG:Debug>,build,build-production>
		WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})

	add_custom_target(webpack DEPENDS ${PROJECT_SOURCE_DIR}/node_modules
		${webpack_output})
endif()

set(CONFIG_PATH "${CMAKE_INSTALL_FULL_SYSCONFDIR}")

set(ALPHAFILL_DATA_DIR
	"${CMAKE_INSTALL_FULL_DATADIR}/alphafill"
	CACHE PATH "The directory where AlphaFill data files are to be stored")

configure_file(${PROJECT_SOURCE_DIR}/src/config.hpp.in
	${PROJECT_SOURCE_DIR}/src/config.hpp)
configure_file(${PROJECT_SOURCE_DIR}/alphafill.conf-dist.in alphafill.conf-dist
	@ONLY)

# pop the BUILD_TESTING flag
set(BUILD_TESTING ${SAVE_BUILD_TESTING})

list(
	APPEND
	SOURCES
	${PROJECT_SOURCE_DIR}/src/config.hpp
	${PROJECT_SOURCE_DIR}/src/alphafill.cpp
	${PROJECT_SOURCE_DIR}/src/blast.cpp
	${PROJECT_SOURCE_DIR}/src/ligands.cpp
	${PROJECT_SOURCE_DIR}/src/matrix.cpp
	${PROJECT_SOURCE_DIR}/src/main.cpp
	${PROJECT_SOURCE_DIR}/src/utilities.cpp
	${PROJECT_SOURCE_DIR}/src/validate.cpp)

if(BUILD_WEB_APPLICATION)
	list(
		APPEND
		SOURCES
		${PROJECT_SOURCE_DIR}/src/bsd-closefrom.c
		${PROJECT_SOURCE_DIR}/src/db-connection.cpp
		${PROJECT_SOURCE_DIR}/src/data-service.cpp
		${PROJECT_SOURCE_DIR}/src/server.cpp
		${PROJECT_SOURCE_DIR}/src/structure.cpp
		${PROJECT_SOURCE_DIR}/src/https-client.cpp)
endif()

# The alphafill executable
add_executable(alphafill ${SOURCES})

find_file(HAVE_DIRENT_H dirent.h)

if(HAVE_DIRENT_H)
	target_compile_definitions(alphafill PUBLIC HAVE_DIRENT_H)
endif()

find_file(HAVE_FCNTL_H fcntl.h)

if(HAVE_FCNTL_H)
	target_compile_definitions(alphafill PUBLIC HAVE_FCNTL_H)
endif()

check_function_exists(dirfd HAVE_DIRFD)

if(HAVE_DIRFD)
	target_compile_definitions(alphafill PUBLIC HAVE_DIRFD)
endif()

check_symbol_exists(closefrom "unistd.h" HAVE_CLOSEFROM)

if(HAVE_CLOSEFROM)
	target_compile_definitions(alphafill PUBLIC HAVE_CLOSEFROM)
endif()

check_symbol_exists(sysconf "unistd.h" HAVE_SYSCONF)

if(HAVE_SYSCONF)
	target_compile_definitions(alphafill PUBLIC HAVE_SYSCONF)
endif()

target_link_libraries(alphafill libmcfp::libmcfp cifpp::cifpp zeep::zeep)

if(BUILD_WEB_APPLICATION)
	target_compile_definitions(alphafill PUBLIC BUILD_WEB_APPLICATION)
	add_dependencies(alphafill webpack)

	target_link_libraries(alphafill libpqxx::pqxx OpenSSL::SSL OpenSSL::Crypto)
endif()

target_compile_definitions(alphafill
	PUBLIC ALPHAFILL_DATA_DIR="${ALPHAFILL_DATA_DIR}")

if(USE_RSRC)
	message("Using resources compiled with ${MRC_EXECUTABLE}")
	add_compile_definitions(USE_RSRC WEBAPP_USES_RESOURCES)

	list(APPEND RESOURCES ${PROJECT_SOURCE_DIR}/rsrc/mmcif_af.dic
		${PROJECT_SOURCE_DIR}/rsrc/mmcif_ma.dic)

	if(CIFPP_SHARE_DIR AND EXISTS "${CIFPP_SHARE_DIR}/mmcif_pdbx.dic")
		list(APPEND RESOURCES ${CIFPP_SHARE_DIR}/mmcif_pdbx.dic)
	endif()

	if(BUILD_WEB_APPLICATION)
		list(APPEND RESOURCES ${PROJECT_SOURCE_DIR}/docroot/
			${PROJECT_SOURCE_DIR}/db-schema.sql
			${PROJECT_SOURCE_DIR}/scripts/refine.mcr)
	endif()

	mrc_target_resources(alphafill ${RESOURCES})
endif()

# installation
install(TARGETS alphafill RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

install(FILES ${PROJECT_SOURCE_DIR}/af-ligands.cif
	DESTINATION ${ALPHAFILL_DATA_DIR})

install(
	CODE "
	if(EXISTS \"${CMAKE_INSTALL_FULL_SYSCONFDIR}/alphafill.conf\")
		file(INSTALL ${CMAKE_BINARY_DIR}/alphafill.conf-dist
		DESTINATION ${CMAKE_INSTALL_FULL_SYSCONFDIR})
	else()
		file(INSTALL ${CMAKE_BINARY_DIR}/alphafill.conf-dist
		DESTINATION ${CMAKE_INSTALL_FULL_SYSCONFDIR}
		RENAME alphafill.conf)
	endif()

	message(NOTICE \"Please edit the ${CMAKE_INSTALL_FULL_SYSCONFDIR}/alphafill.conf file.\")
")

# If we do not use resources, install the resource data in a share directory
if(NOT USE_RSRC)
	install(
		FILES ${CIFPP_SHARE_DIR}/mmcif_pdbx.dic
		${PROJECT_SOURCE_DIR}/rsrc/mmcif_af.dic
		${PROJECT_SOURCE_DIR}/rsrc/mmcif_ma.dic
		DESTINATION ${ALPHAFILL_DATA_DIR})

	if(BUILD_WEB_APPLICATION)
		install(FILES ${PROJECT_SOURCE_DIR}/db-schema.sql
			${PROJECT_SOURCE_DIR}/scripts/refine.mcr
			DESTINATION ${ALPHAFILL_DATA_DIR})

		install(DIRECTORY ${PROJECT_SOURCE_DIR}/docroot
			DESTINATION ${ALPHAFILL_DATA_DIR})
	endif()
endif()

if(BUILD_DOCUMENTATION)
	add_subdirectory(docs)
endif()

if(BUILD_TESTING)
	# simply run alphafill jobs on our test data start by generating a new config
	# file
	configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/alphafill-test.conf.in
		${CMAKE_CURRENT_BINARY_DIR}/alphafill-test.conf @ONLY)

	# First create an index
	add_test(NAME create_index_test COMMAND alphafill create-index
		--config=alphafill-test.conf)

	add_test(
		NAME run_process_test
		COMMAND
		alphafill process --config=alphafill-test.conf
		${CMAKE_CURRENT_SOURCE_DIR}/test/afdb-v4/P2/AF-P29373-F1-model_v4.cif.gz
		AF-P29373-F1-filled.cif.gz
		WORKING_DIRECTORY)

	add_executable(validate-test-files
		${PROJECT_SOURCE_DIR}/test/validate-test-files.cpp)
	target_link_libraries(validate-test-files libmcfp::libmcfp cifpp::cifpp
		zeep::zeep)

	add_test(NAME validate_files_test
		COMMAND validate-test-files AF-P29373-F1-filled.cif.gz
		AF-P29373-F1-filled.cif.json --pdb-id=1CBS --pdb-asym-id=A)

	# If cifpp was not found, it means it was not installed so we might not have a
	# proper CIFPP_SHARE_DIR and the the tests will fail
	if(NOT cifpp_FOUND)
		set_property(TEST run_process_test
			PROPERTY ENVIRONMENT "LIBCIFPP_DATA_DIR=${CIFPP_SHARE_DIR}")
		set_property(TEST validate_files_test
			PROPERTY ENVIRONMENT "LIBCIFPP_DATA_DIR=${CIFPP_SHARE_DIR}")
	endif()
endif()
