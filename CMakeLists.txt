# SPDX-License-Identifier: BSD-2-Clause

# Copyright (c) 2021 NKI/AVL, Netherlands Cancer Institute

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# 1. Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimer
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.

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

cmake_minimum_required(VERSION 3.21)

# set the project name
project(alphafill VERSION 2.1.0 LANGUAGES C CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(GNUInstallDirs)
include(CheckFunctionExists)
include(CheckIncludeFiles)
include(CheckLibraryExists)
include(CMakePackageConfigHelpers)
include(GenerateExportHeader)
include(CTest)
include(FetchContent)

set(CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(BUILD_DOCUMENTATION "Build the documentation" OFF)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-unused-parameter -Wno-missing-field-initializers")
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
endif()

# Create a revision file, containing the current git version info
include(VersionString)
write_version_header("${PROJECT_SOURCE_DIR}/src/")

find_package(zeep 6.0.9 QUIET)
if (NOT zeep_FOUND)
	FetchContent_Declare(zeep GIT_REPOSITORY https://github.com/mhekkel/libzeep.git GIT_TAG v6.0.9)
	FetchContent_MakeAvailable(zeep)
endif()

find_package(cifpp 6 QUIET)
if (NOT cifpp_FOUND)
	FetchContent_Declare(cifpp GIT_REPOSITORY https://github.com/PDB-REDO/libcifpp.git GIT_TAG 2692f2c)
	FetchContent_MakeAvailable(cifpp)
endif()

find_package(libmcfp 1.2.4 QUIET)
if (NOT libmcfp_FOUND)
	FetchContent_Declare(libmcfp GIT_REPOSITORY https://github.com/mhekkel/libcifpp.git GIT_TAG v1.2.4)
	FetchContent_MakeAvailable(libmcfp)
endif()

find_package(libpqxx 7.8.0 QUIET)
if (NOT libpqxx_FOUND)
	FetchContent_Declare(libpqxx GIT_REPOSITORY https://github.com/jtv/libpqxx GIT_TAG 7.8.0)
	FetchContent_MakeAvailable(libpqxx)
endif()

find_package(OpenSSL REQUIRED)

find_program(YARN yarn REQUIRED)

# yarn rules for javascripts
set(webpack_input
	${PROJECT_SOURCE_DIR}/webapp/index.js
	${PROJECT_SOURCE_DIR}/webapp/model.js
	${PROJECT_SOURCE_DIR}/webapp/optimized.js
	${PROJECT_SOURCE_DIR}/webapp/lists.js
	${PROJECT_SOURCE_DIR}/webapp/wait.js
	${PROJECT_SOURCE_DIR}/webapp/molstar.tsx
	${PROJECT_SOURCE_DIR}/scss/pdb-redo-bootstrap.scss
)

set(webpack_output ${PROJECT_SOURCE_DIR}/docroot/scripts/index.js)

add_custom_command(
	OUTPUT ${webpack_output}
	BYPRODUCTS
	DEPENDS ${webpack_input}
	COMMAND ${YARN} run $<IF:$<CONFIG:Debug>,build,build-production>
	WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})

add_custom_target(webpack DEPENDS ${webpack_output})

find_file(HAVE_DIRENT_H dirent.h)
find_file(HAVE_FCNTL_H fcntl.h)
check_function_exists(dirfd HAVE_DIRFD)
find_file(HAVE_PROC_PID '/proc/self/fd/0')
check_function_exists(sysconf HAVE_SYSCONF)

# The alphafill executable
add_executable(alphafill
	${PROJECT_SOURCE_DIR}/src/alphafill.cpp
	${PROJECT_SOURCE_DIR}/src/blast.cpp
	${PROJECT_SOURCE_DIR}/src/bsd-closefrom.c
	${PROJECT_SOURCE_DIR}/src/ligands.cpp
	${PROJECT_SOURCE_DIR}/src/matrix.cpp
	${PROJECT_SOURCE_DIR}/src/main.cpp
	${PROJECT_SOURCE_DIR}/src/utilities.cpp

	${PROJECT_SOURCE_DIR}/src/db-connection.cpp
	${PROJECT_SOURCE_DIR}/src/data-service.cpp
	${PROJECT_SOURCE_DIR}/src/server.cpp
	${PROJECT_SOURCE_DIR}/src/structure.cpp
	${PROJECT_SOURCE_DIR}/src/validate.cpp

	${PROJECT_SOURCE_DIR}/src/https-client.cpp
)

target_compile_definitions(alphafill PUBLIC
	HAVE_DIRENT_H=${HAVE_DIRENT_H}
	HAVE_FCNTL_H=${HAVE_FCNTL_H}
	HAVE_DIRFD=${HAVE_DIRFD}
	HAVE_PROC_PID=${HAVE_PROC_PID}
	HAVE_SYSCONF=${HAVE_SYSCONF}
)

add_dependencies(alphafill webpack)

target_link_libraries(alphafill libmcfp::libmcfp cifpp::cifpp zeep::zeep libpqxx::pqxx OpenSSL::SSL OpenSSL::Crypto)

# Optionally use mrc to create resources
find_package(Mrc QUIET)

if(MRC_FOUND)
	option(USE_RSRC "Use mrc to create resources" ON)
else()
	message(STATUS "Not using resources since mrc was not found")
endif()

if(USE_RSRC)
	message("Using resources compiled with ${MRC_EXECUTABLE}")
	add_compile_definitions(USE_RSRC WEBAPP_USES_RESOURCES)

	mrc_target_resources(alphafill

		${PROJECT_SOURCE_DIR}/docroot/
		${PROJECT_SOURCE_DIR}/db-schema.sql
		${PROJECT_SOURCE_DIR}/scripts/refine.mcr

		${CIFPP_SHARE_DIR}/mmcif_pdbx.dic
		${PROJECT_SOURCE_DIR}/rsrc/mmcif_af.dic
		${PROJECT_SOURCE_DIR}/rsrc/mmcif_ma.dic)
else()
	set(ALPHAFILL_DATA_DIR $<IF:$<CONFIG:Debug>,${PROJECT_SOURCE_DIR}/rsrc,${CMAKE_INSTALL_FULL_DATADIR}/alphafill>)
	add_compile_definitions(ALPHAFILL_DATA_DIR="${ALPHAFILL_DATA_DIR}")
endif()

# installation
install(TARGETS alphafill RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

list(APPEND hardlinks process create-index rebuild-db server)

foreach(hardlink ${hardlinks})
	install(CODE "set(link_from \"${CMAKE_INSTALL_FULL_BINDIR}/alphafill\")")
	install(CODE "set(link_to \"${CMAKE_INSTALL_FULL_BINDIR}/alphafill-${hardlink}\")")
	install(CODE [[
message(STATUS "Creating hard link from ${link_from} to ${link_to}")
execute_process(COMMAND ${CMAKE_COMMAND} -E create_hardlink ${link_from} ${link_to})
]])
endforeach()

# If we do not use resources, install the resource data in a share directory
if(NOT USE_RSRC)
	install(FILES
		${PROJECT_SOURCE_DIR}/db-schema.sql
		${PROJECT_SOURCE_DIR}/scripts/refine.mcr

		${CIFPP_SHARE_DIR}/mmcif_pdbx.dic
		${PROJECT_SOURCE_DIR}/rsrc/mmcif_af.dic
		${PROJECT_SOURCE_DIR}/rsrc/mmcif_ma.dic

		DESTINATION ${ALPHAFILL_DATA_DIR}
	)

	install(DIRECTORY
		${PROJECT_SOURCE_DIR}/docroot
		DESTINATION ${ALPHAFILL_DATA_DIR}
	)
endif()

if(BUILD_DOCUMENTATION)
	add_subdirectory(docs)
endif()