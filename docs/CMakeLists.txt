find_package(Sphinx REQUIRED)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in ${CMAKE_CURRENT_SOURCE_DIR}/conf.py @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/man/alphafill.conf.rst.in ${CMAKE_CURRENT_SOURCE_DIR}/man/alphafill.conf.rst @ONLY)

set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR})
# set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/sphinx)
set(SPHINX_BUILD ${PROJECT_SOURCE_DIR}/docroot)

add_custom_target("Sphinx-${PROJECT_NAME}-html" ALL
	COMMAND ${SPHINX_EXECUTABLE} -b dirhtml
	${SPHINX_SOURCE} ${SPHINX_BUILD}
	COMMENT "Generating documentation with Sphinx")

add_custom_target("Sphinx-${PROJECT_NAME}-man" ALL
	COMMAND ${SPHINX_EXECUTABLE} -b man
	${SPHINX_SOURCE} ${CMAKE_CURRENT_BINARY_DIR}/man
	COMMENT "Generating manual pages with Sphinx")

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/man/alphafill.1
	${CMAKE_CURRENT_BINARY_DIR}/man/alphafill-api.1
	${CMAKE_CURRENT_BINARY_DIR}/man/alphafill-create-index.1
	${CMAKE_CURRENT_BINARY_DIR}/man/alphafill-process.1
	${CMAKE_CURRENT_BINARY_DIR}/man/alphafill-rebuild-db.1
	${CMAKE_CURRENT_BINARY_DIR}/man/alphafill-server.1
	DESTINATION ${CMAKE_INSTALL_DATADIR}/man/man1)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/man/alphafill.conf.5
	DESTINATION ${CMAKE_INSTALL_DATADIR}/man/man5)
