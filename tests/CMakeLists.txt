add_executable(test_bsm test_bsm.c)
target_link_libraries(test_bsm PRIVATE WonderVol::wondervol)
add_test(NAME bsm_tests COMMAND test_bsm)

add_executable(test_iv test_iv.c)
target_link_libraries(test_iv PRIVATE WonderVol::wondervol)
add_test(NAME iv_tests COMMAND test_iv)

add_executable(test_dispersion test_dispersion.c)
target_link_libraries(test_dispersion PRIVATE WonderVol::wondervol)
add_test(NAME dispersion_tests COMMAND test_dispersion)

add_executable(test_platform test_platform.c)
target_link_libraries(test_platform PRIVATE WonderVol::wondervol)
add_test(NAME platform_tests COMMAND test_platform)
