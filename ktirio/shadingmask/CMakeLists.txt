add_custom_target(pyshadingmask)
feelpp_add_pymodule(shadingmask SRCS shadingmask.cpp DESTINATION shadingmask )
target_include_directories(_shadingmask PRIVATE ${CMAKE_SOURCE_DIR}/src)
add_dependencies(pyshadingmask _shadingmask)

install(FILES  __init__.py DESTINATION ${FEELPP_PYTHON_MODULE_PATH}/shadingmask)