add_executable(bench_iv bench_iv.c)
target_link_libraries(bench_iv PRIVATE Dispersion::Scanner)

add_executable(bench_dispersion bench_dispersion.c)
target_link_libraries(bench_dispersion PRIVATE Dispersion::Scanner)

add_executable(bench_neon_scalar bench_neon_scalar.c)
target_link_libraries(bench_neon_scalar PRIVATE Dispersion::Scanner)