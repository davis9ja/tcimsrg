cmake_minimum_required(VERSION 3.11)

# set( CMAKE_C_COMPILER "icpc" )
set( CMAKE_CXX_COMPILER "g++" )
# set( MPI_C_COMPILER "mpiicpc" )
# set( MPI_CXX_COMPILER "mpiicpc" )
set(PROJECT_DIR .)
set(MKLROOT $ENV{MKLROOT})

# set the project name
project("imsrg" LANGUAGES CXX)

#find_package(taco CONFIG REQUIRED)
find_package(OpenMP REQUIRED)
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# add the executable
include_directories(/mnt/home/daviso53/taco/include ${PROJECT_DIR}/include /mnt/home/daviso53/Research/dmcpp/ ${MKLROOT}/include)
file(GLOB SOURCES ${PROJECT_DIR}/src/*.cpp)
add_executable(solve_imsrg_test ${SOURCES})
target_link_libraries(solve_imsrg_test /mnt/home/daviso53/taco/build/lib/libtaco.so /mnt/home/daviso53/Research/dmcpp/libdmc.a -lstdc++ -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_openmpi_ilp64.a -Wl,--end-group  -fopenmp -lpthread -lm -ldl)
add_compile_options(-DMKL_ILP64  -m64 -O1)

