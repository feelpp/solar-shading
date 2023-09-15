set(PYFILES conftest.py test_surface.py test_volume.py)
install(FILES ${PYFILES} DESTINATION ${FEELPP_PYTHON_MODULE_PATH}/shadingmask)

foreach( test  surface volume )
add_test (NAME feelpp_sm_python-${test}
  COMMAND ${CMAKE_COMMAND} -E env
        LD_LIBRARY_PATH=:${CMAKE_INSTALL_PREFIX}/lib:$ENV{LD_LIBRARY_PATH}
        PYTHONPATH=${CMAKE_INSTALL_PREFIX}/${FEELPP_PYTHON_MODULE_PATH}:$ENV{PYTHONPATH} #${CMAKE_CURRENT_BINARY_DIR}/..:${CMAKE_CURRENT_SOURCE_DIR}/..
        ${PYTHON_EXECUTABLE} -m pytest -k ${test} -s -vv ${CMAKE_CURRENT_SOURCE_DIR}
)
endforeach()