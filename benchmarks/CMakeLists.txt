add_executable(bench_iv bench_iv.c)
target_link_libraries(bench_iv PRIVATE WonderVol::wondervol)

add_executable(bench_dispersion bench_dispersion.c)
target_link_libraries(bench_dispersion PRIVATE WonderVol::wondervol)

add_executable(bench_simd_scalar bench_neon_scalar.c)
target_link_libraries(bench_simd_scalar PRIVATE WonderVol::wondervol)
