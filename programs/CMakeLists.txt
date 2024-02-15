if (WIN32)
    set(ONEAPI_ROOT "C:/Program\ Files\ (x86)/Intel/oneAPI")
endif()
if(UNIX)
    set(MUMPS_ROOT "/mnt/c/Users/Compucter/libs/mumps-main")
    set(OPENBLAS_ROOT "/mnt/c/Users/Compucter/libs/OpenBLAS-develop")
endif()

set(COMMON_EXECUTABLE data_formats.c data_formats.h mmio.c mmio.h)

# These programs use internal metis data structures.
include_directories(../libmetis)
link_directories(/home/karypis/local/lib)
# Build program.
add_executable(demo_metis demo_metis.c ${COMMON_EXECUTABLE})
if (WIN32)
    add_executable(demo_metis_pardiso demo_metis_pardiso.c ${COMMON_EXECUTABLE})

    target_include_directories(demo_metis_pardiso 
        PUBLIC ${ONEAPI_ROOT}/mkl/latest/include
    )

    target_link_libraries(demo_metis_pardiso
        PUBLIC ${ONEAPI_ROOT}/mkl/latest/lib/intel64/mkl_intel_lp64.lib
        PUBLIC ${ONEAPI_ROOT}/mkl/latest/lib/intel64/mkl_intel_thread.lib
        PUBLIC ${ONEAPI_ROOT}/mkl/latest/lib/intel64/mkl_core.lib
        PUBLIC ${ONEAPI_ROOT}/compiler/latest/windows/compiler/lib/intel64_win/libiomp5md.lib
    )
endif()
if(UNIX)
    add_executable(demo_metis_mumps demo_metis_mumps.c ${COMMON_EXECUTABLE})

    target_include_directories(demo_metis_mumps
        PUBLIC ${MUMPS_ROOT}/cache/MUMPS-5.3.5-src/include
        PUBLIC ${OPENBLAS_ROOT}/include
        )

    target_link_libraries(demo_metis_mumps
        PUBLIC ${MUMPS_ROOT}/build/libdmumps.a
        PUBLIC ${MUMPS_ROOT}/build/libmumps_common.a
        PUBLIC ${MUMPS_ROOT}/build/libpord.a
        PUBLIC ${MUMPS_ROOT}/build/lib/libblas.a
        PUBLIC ${MUMPS_ROOT}/build/lib/liblapack.a
        PUBLIC ${MUMPS_ROOT}/build/_deps/mumps-build/libseq/libmpiseq.a
        PUBLIC ${OPENBLAS_ROOT}/libopenblas_zenp-r0.3.19.dev.a
        -lm -lgfortran -lpthread
        )
endif()

foreach(prog demo_metis demo_metis_mumps)
  target_link_libraries(${prog} metis)
#  target_link_libraries(${prog} metis profiler)
endforeach(prog)

if(METIS_INSTALL)
  install(TARGETS demo_metis demo_metis_mumps
    RUNTIME DESTINATION bin)
endif()

# Try to find subversion revision.
set(SVNREV "")
file(TO_NATIVE_PATH ${PROJECT_SOURCE_DIR}/.svn svn_dir)
if(IS_DIRECTORY ${svn_dir})
  include(FindSubversion)
  if(Subversion_FOUND)
    Subversion_WC_INFO(${PROJECT_SOURCE_DIR} metis)
    set(SVNREV ${metis_WC_REVISION})
  endif(Subversion_FOUND)
endif()
add_definitions(-DSVNINFO="${SVNREV}")
