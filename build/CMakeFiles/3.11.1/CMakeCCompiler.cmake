set(CMAKE_C_COMPILER "/opt/software/GCCcore/6.4.0/bin/cc")
set(CMAKE_C_COMPILER_ARG1 "")
set(CMAKE_C_COMPILER_ID "GNU")
set(CMAKE_C_COMPILER_VERSION "6.4.0")
set(CMAKE_C_COMPILER_VERSION_INTERNAL "")
set(CMAKE_C_COMPILER_WRAPPER "")
set(CMAKE_C_STANDARD_COMPUTED_DEFAULT "11")
set(CMAKE_C_COMPILE_FEATURES "c_std_90;c_function_prototypes;c_std_99;c_restrict;c_variadic_macros;c_std_11;c_static_assert")
set(CMAKE_C90_COMPILE_FEATURES "c_std_90;c_function_prototypes")
set(CMAKE_C99_COMPILE_FEATURES "c_std_99;c_restrict;c_variadic_macros")
set(CMAKE_C11_COMPILE_FEATURES "c_std_11;c_static_assert")

set(CMAKE_C_PLATFORM_ID "Linux")
set(CMAKE_C_SIMULATE_ID "")
set(CMAKE_C_SIMULATE_VERSION "")



set(CMAKE_AR "/opt/software/binutils/2.28-GCCcore-6.4.0/bin/ar")
set(CMAKE_C_COMPILER_AR "/opt/software/GCCcore/6.4.0/bin/gcc-ar")
set(CMAKE_RANLIB "/opt/software/binutils/2.28-GCCcore-6.4.0/bin/ranlib")
set(CMAKE_C_COMPILER_RANLIB "/opt/software/GCCcore/6.4.0/bin/gcc-ranlib")
set(CMAKE_LINKER "/opt/software/binutils/2.28-GCCcore-6.4.0/bin/ld")
set(CMAKE_COMPILER_IS_GNUCC 1)
set(CMAKE_C_COMPILER_LOADED 1)
set(CMAKE_C_COMPILER_WORKS TRUE)
set(CMAKE_C_ABI_COMPILED TRUE)
set(CMAKE_COMPILER_IS_MINGW )
set(CMAKE_COMPILER_IS_CYGWIN )
if(CMAKE_COMPILER_IS_CYGWIN)
  set(CYGWIN 1)
  set(UNIX 1)
endif()

set(CMAKE_C_COMPILER_ENV_VAR "CC")

if(CMAKE_COMPILER_IS_MINGW)
  set(MINGW 1)
endif()
set(CMAKE_C_COMPILER_ID_RUN 1)
set(CMAKE_C_SOURCE_FILE_EXTENSIONS c;m)
set(CMAKE_C_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_C_LINKER_PREFERENCE 10)

# Save compiler ABI information.
set(CMAKE_C_SIZEOF_DATA_PTR "8")
set(CMAKE_C_COMPILER_ABI "ELF")
set(CMAKE_C_LIBRARY_ARCHITECTURE "")

if(CMAKE_C_SIZEOF_DATA_PTR)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_C_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_C_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_C_COMPILER_ABI}")
endif()

if(CMAKE_C_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()

set(CMAKE_C_CL_SHOWINCLUDES_PREFIX "")
if(CMAKE_C_CL_SHOWINCLUDES_PREFIX)
  set(CMAKE_CL_SHOWINCLUDES_PREFIX "${CMAKE_C_CL_SHOWINCLUDES_PREFIX}")
endif()





set(CMAKE_C_IMPLICIT_LINK_LIBRARIES "gcc;gcc_s;c;gcc;gcc_s")
set(CMAKE_C_IMPLICIT_LINK_DIRECTORIES "/opt/software/cuDNN/7.5.0.56-CUDA-10.0.130/lib64;/opt/software/CUDA/10.0.130/lib64;/opt/software/libffi/3.2.1-GCCcore-6.4.0/lib64;/opt/software/GCCcore/6.4.0/lib64;/opt/software/GCCcore/6.4.0/lib/gcc/x86_64-pc-linux-gnu/6.4.0;/lib64;/usr/lib64;/opt/software/Python/3.6.4-foss-2018a/lib;/opt/software/libffi/3.2.1-GCCcore-6.4.0/lib;/opt/software/GMP/6.1.2-GCCcore-6.4.0/lib;/opt/software/SQLite/3.21.0-GCCcore-6.4.0/lib;/opt/software/libreadline/7.0-GCCcore-6.4.0/lib;/opt/software/ncurses/6.0-GCCcore-6.4.0/lib;/opt/software/imkl/2018.1.163-gompi-2018a/mkl/lib/intel64;/opt/software/imkl/2018.1.163-gompi-2018a/lib/intel64;/opt/software/Boost/1.67.0-foss-2018a/lib;/opt/software/zlib/1.2.11-GCCcore-6.4.0/lib;/opt/software/bzip2/1.0.6-GCCcore-6.4.0/lib;/opt/software/ScaLAPACK/2.0.2-gompi-2018a-OpenBLAS-0.2.20/lib;/opt/software/FFTW/3.3.7-gompi-2018a/lib;/opt/software/OpenBLAS/0.2.20-GCC-6.4.0-2.28/lib;/opt/software/OpenMPI/2.1.2-GCC-6.4.0-2.28/lib;/opt/software/binutils/2.28-GCCcore-6.4.0/lib;/opt/software/GCCcore/6.4.0/lib;/opt/software/Java/1.8.0_152/lib;/opt/software/Tcl/8.6.8-GCCcore-6.4.0/lib;/opt/software/tbb/2018_U3-GCCcore-6.4.0/build/linux_intel64_gcc_cc6.4.0_libc2.17_kernel3.10.0_release")
set(CMAKE_C_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
