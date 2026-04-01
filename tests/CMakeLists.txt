add_executable(test_bsm test_bsm.c)
target_link_libraries(test_bsm PRIVATE Dispersion::Scanner)
add_test(NAME bsm_tests COMMAND test_bsm)

add_executable(test_iv test_iv.c)
target_link_libraries(test_iv PRIVATE Dispersion::Scanner)
add_test(NAME iv_tests COMMAND test_iv)

add_executable(test_dispersion test_dispersion.c)
target_link_libraries(test_dispersion PRIVATE Dispersion::Scanner)
add_test(NAME dispersion_tests COMMAND test_dispersion)