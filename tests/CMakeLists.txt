enable_testing()

find_package(Python3 REQUIRED COMPONENTS Interpreter)

add_executable(test_hydrophone_unpacker
    test_hydrophone_unpacker.cpp
    test_utils.cpp
)
target_link_libraries(test_hydrophone_unpacker PRIVATE sonar_core)
add_test(NAME test_hydrophone_unpacker COMMAND test_hydrophone_unpacker)

add_executable(test_cma_equalizer
    test_cma_equalizer.cpp
    test_utils.cpp
)
target_link_libraries(test_cma_equalizer PRIVATE sonar_core)
add_test(NAME test_cma_equalizer COMMAND test_cma_equalizer)

add_executable(test_doppler_spectrum
    test_doppler_spectrum.cpp
    test_utils.cpp
)
target_link_libraries(test_doppler_spectrum PRIVATE sonar_core)
add_test(NAME test_doppler_spectrum COMMAND test_doppler_spectrum)

add_executable(generate_test_data
    generate_test_data.cpp
    test_utils.cpp
)
target_link_libraries(generate_test_data PRIVATE sonar_core)

if(CMAKE_CONFIGURATION_TYPES)
    set(PYTHON_MODULE_DIR ${CMAKE_BINARY_DIR}/python/$<CONFIG>)
else()
    set(PYTHON_MODULE_DIR ${CMAKE_BINARY_DIR}/python)
endif()

add_test(
    NAME test_python_integration
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_python.py
    WORKING_DIRECTORY ${PYTHON_MODULE_DIR}
)

set_tests_properties(test_python_integration PROPERTIES
    DEPENDS sonar_py
    ENVIRONMENT "PYTHONPATH=${PYTHON_MODULE_DIR}"
)
