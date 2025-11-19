find_package(CUDAToolkit REQUIRED)
enable_language(CUDA)

add_library(forrf_lib
    fftpack.f90
    sacio.f90
    deconit.f90
    deconit.cu
)

target_link_libraries(forrf_lib PRIVATE CUDA::cufft)
target_compile_options(forrf_lib PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--extended-lambda>)




add_executable(forrf_test main.f90)
target_link_libraries(forrf_test PRIVATE forrf_lib)

file(COPY ${CMAKE_SOURCE_DIR}/test_R.sac ${CMAKE_SOURCE_DIR}/test_T.sac ${CMAKE_SOURCE_DIR}/test_Z.sac DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME RunDeconTest COMMAND forrf_test)
set_tests_properties(RunDeconTest PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})


