cmake_minimum_required(VERSION 3.11)

# set( CMAKE_C_COMPILER "icpc" )
set( CMAKE_CXX_COMPILER "g++" )
# set( MPI_C_COMPILER "mpiicpc" )
# set( MPI_CXX_COMPILER "mpiicpc" )
set(PROJECT_DIR .)

# set the project name
project("imsrg" LANGUAGES CXX)

#find_package(taco CONFIG REQUIRED)
find_package(OpenMP REQUIRED)
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()


# add the executable
include_directories(/mnt/home/daviso53/taco/include ${PROJECT_DIR}/include /mnt/home/daviso53/Research/dmcpp/)
file(GLOB SOURCES ${PROJECT_DIR}/src/*.cpp)
add_executable(solve_imsrg ${SOURCES})
target_link_libraries(solve_imsrg /mnt/home/daviso53/taco/build/lib/libtaco.so /mnt/home/daviso53/Research/dmcpp/libdmc.a)
