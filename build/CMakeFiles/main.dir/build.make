# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/software/CMake/3.20.1-GCCcore-10.3.0/bin/cmake

# The command to remove a file.
RM = /opt/software/CMake/3.20.1-GCCcore-10.3.0/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /mnt/home/daviso53/Research/tcimsrg

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /mnt/home/daviso53/Research/tcimsrg/build

# Include any dependencies generated for this target.
include CMakeFiles/main.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/main.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/main.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/main.dir/flags.make

CMakeFiles/main.dir/src/hwrapper.cpp.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/src/hwrapper.cpp.o: ../src/hwrapper.cpp
CMakeFiles/main.dir/src/hwrapper.cpp.o: CMakeFiles/main.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/home/daviso53/Research/tcimsrg/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/main.dir/src/hwrapper.cpp.o"
	/opt/software/GCCcore/10.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/main.dir/src/hwrapper.cpp.o -MF CMakeFiles/main.dir/src/hwrapper.cpp.o.d -o CMakeFiles/main.dir/src/hwrapper.cpp.o -c /mnt/home/daviso53/Research/tcimsrg/src/hwrapper.cpp

CMakeFiles/main.dir/src/hwrapper.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main.dir/src/hwrapper.cpp.i"
	/opt/software/GCCcore/10.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/home/daviso53/Research/tcimsrg/src/hwrapper.cpp > CMakeFiles/main.dir/src/hwrapper.cpp.i

CMakeFiles/main.dir/src/hwrapper.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main.dir/src/hwrapper.cpp.s"
	/opt/software/GCCcore/10.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/home/daviso53/Research/tcimsrg/src/hwrapper.cpp -o CMakeFiles/main.dir/src/hwrapper.cpp.s

CMakeFiles/main.dir/src/main.cpp.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/src/main.cpp.o: ../src/main.cpp
CMakeFiles/main.dir/src/main.cpp.o: CMakeFiles/main.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/home/daviso53/Research/tcimsrg/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/main.dir/src/main.cpp.o"
	/opt/software/GCCcore/10.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/main.dir/src/main.cpp.o -MF CMakeFiles/main.dir/src/main.cpp.o.d -o CMakeFiles/main.dir/src/main.cpp.o -c /mnt/home/daviso53/Research/tcimsrg/src/main.cpp

CMakeFiles/main.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main.dir/src/main.cpp.i"
	/opt/software/GCCcore/10.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/home/daviso53/Research/tcimsrg/src/main.cpp > CMakeFiles/main.dir/src/main.cpp.i

CMakeFiles/main.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main.dir/src/main.cpp.s"
	/opt/software/GCCcore/10.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/home/daviso53/Research/tcimsrg/src/main.cpp -o CMakeFiles/main.dir/src/main.cpp.s

CMakeFiles/main.dir/src/occupation_factors.cpp.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/src/occupation_factors.cpp.o: ../src/occupation_factors.cpp
CMakeFiles/main.dir/src/occupation_factors.cpp.o: CMakeFiles/main.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/home/daviso53/Research/tcimsrg/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/main.dir/src/occupation_factors.cpp.o"
	/opt/software/GCCcore/10.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/main.dir/src/occupation_factors.cpp.o -MF CMakeFiles/main.dir/src/occupation_factors.cpp.o.d -o CMakeFiles/main.dir/src/occupation_factors.cpp.o -c /mnt/home/daviso53/Research/tcimsrg/src/occupation_factors.cpp

CMakeFiles/main.dir/src/occupation_factors.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main.dir/src/occupation_factors.cpp.i"
	/opt/software/GCCcore/10.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/home/daviso53/Research/tcimsrg/src/occupation_factors.cpp > CMakeFiles/main.dir/src/occupation_factors.cpp.i

CMakeFiles/main.dir/src/occupation_factors.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main.dir/src/occupation_factors.cpp.s"
	/opt/software/GCCcore/10.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/home/daviso53/Research/tcimsrg/src/occupation_factors.cpp -o CMakeFiles/main.dir/src/occupation_factors.cpp.s

CMakeFiles/main.dir/src/pairinghamiltonian.cpp.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/src/pairinghamiltonian.cpp.o: ../src/pairinghamiltonian.cpp
CMakeFiles/main.dir/src/pairinghamiltonian.cpp.o: CMakeFiles/main.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/home/daviso53/Research/tcimsrg/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/main.dir/src/pairinghamiltonian.cpp.o"
	/opt/software/GCCcore/10.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/main.dir/src/pairinghamiltonian.cpp.o -MF CMakeFiles/main.dir/src/pairinghamiltonian.cpp.o.d -o CMakeFiles/main.dir/src/pairinghamiltonian.cpp.o -c /mnt/home/daviso53/Research/tcimsrg/src/pairinghamiltonian.cpp

CMakeFiles/main.dir/src/pairinghamiltonian.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main.dir/src/pairinghamiltonian.cpp.i"
	/opt/software/GCCcore/10.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/home/daviso53/Research/tcimsrg/src/pairinghamiltonian.cpp > CMakeFiles/main.dir/src/pairinghamiltonian.cpp.i

CMakeFiles/main.dir/src/pairinghamiltonian.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main.dir/src/pairinghamiltonian.cpp.s"
	/opt/software/GCCcore/10.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/home/daviso53/Research/tcimsrg/src/pairinghamiltonian.cpp -o CMakeFiles/main.dir/src/pairinghamiltonian.cpp.s

CMakeFiles/main.dir/src/white.cpp.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/src/white.cpp.o: ../src/white.cpp
CMakeFiles/main.dir/src/white.cpp.o: CMakeFiles/main.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/home/daviso53/Research/tcimsrg/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/main.dir/src/white.cpp.o"
	/opt/software/GCCcore/10.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/main.dir/src/white.cpp.o -MF CMakeFiles/main.dir/src/white.cpp.o.d -o CMakeFiles/main.dir/src/white.cpp.o -c /mnt/home/daviso53/Research/tcimsrg/src/white.cpp

CMakeFiles/main.dir/src/white.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main.dir/src/white.cpp.i"
	/opt/software/GCCcore/10.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/home/daviso53/Research/tcimsrg/src/white.cpp > CMakeFiles/main.dir/src/white.cpp.i

CMakeFiles/main.dir/src/white.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main.dir/src/white.cpp.s"
	/opt/software/GCCcore/10.3.0/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/home/daviso53/Research/tcimsrg/src/white.cpp -o CMakeFiles/main.dir/src/white.cpp.s

# Object files for target main
main_OBJECTS = \
"CMakeFiles/main.dir/src/hwrapper.cpp.o" \
"CMakeFiles/main.dir/src/main.cpp.o" \
"CMakeFiles/main.dir/src/occupation_factors.cpp.o" \
"CMakeFiles/main.dir/src/pairinghamiltonian.cpp.o" \
"CMakeFiles/main.dir/src/white.cpp.o"

# External object files for target main
main_EXTERNAL_OBJECTS =

main: CMakeFiles/main.dir/src/hwrapper.cpp.o
main: CMakeFiles/main.dir/src/main.cpp.o
main: CMakeFiles/main.dir/src/occupation_factors.cpp.o
main: CMakeFiles/main.dir/src/pairinghamiltonian.cpp.o
main: CMakeFiles/main.dir/src/white.cpp.o
main: CMakeFiles/main.dir/build.make
main: /mnt/home/daviso53/taco/build/lib/libtaco.so
main: CMakeFiles/main.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/mnt/home/daviso53/Research/tcimsrg/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking CXX executable main"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/main.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/main.dir/build: main
.PHONY : CMakeFiles/main.dir/build

CMakeFiles/main.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/main.dir/cmake_clean.cmake
.PHONY : CMakeFiles/main.dir/clean

CMakeFiles/main.dir/depend:
	cd /mnt/home/daviso53/Research/tcimsrg/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt/home/daviso53/Research/tcimsrg /mnt/home/daviso53/Research/tcimsrg /mnt/home/daviso53/Research/tcimsrg/build /mnt/home/daviso53/Research/tcimsrg/build /mnt/home/daviso53/Research/tcimsrg/build/CMakeFiles/main.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/main.dir/depend
