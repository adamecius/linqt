# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/santiago/Documents/ICN2/Codes/linqt-3.0.0_beta

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/santiago/Documents/ICN2/Codes/linqt-3.0.0_beta/build

# Include any dependencies generated for this target.
include src/CMakeFiles/inline_compute-kpm-spectralOp-nonOrth.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/CMakeFiles/inline_compute-kpm-spectralOp-nonOrth.dir/compiler_depend.make

# Include the progress variables for this target.
include src/CMakeFiles/inline_compute-kpm-spectralOp-nonOrth.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/inline_compute-kpm-spectralOp-nonOrth.dir/flags.make

src/CMakeFiles/inline_compute-kpm-spectralOp-nonOrth.dir/inline_compute-kpm-spectralOp-nonOrth.o: src/CMakeFiles/inline_compute-kpm-spectralOp-nonOrth.dir/flags.make
src/CMakeFiles/inline_compute-kpm-spectralOp-nonOrth.dir/inline_compute-kpm-spectralOp-nonOrth.o: ../src/inline_compute-kpm-spectralOp-nonOrth.cpp
src/CMakeFiles/inline_compute-kpm-spectralOp-nonOrth.dir/inline_compute-kpm-spectralOp-nonOrth.o: src/CMakeFiles/inline_compute-kpm-spectralOp-nonOrth.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/santiago/Documents/ICN2/Codes/linqt-3.0.0_beta/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/inline_compute-kpm-spectralOp-nonOrth.dir/inline_compute-kpm-spectralOp-nonOrth.o"
	cd /home/santiago/Documents/ICN2/Codes/linqt-3.0.0_beta/build/src && /opt/intel/oneapi/compiler/2024.1/bin/icpx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/CMakeFiles/inline_compute-kpm-spectralOp-nonOrth.dir/inline_compute-kpm-spectralOp-nonOrth.o -MF CMakeFiles/inline_compute-kpm-spectralOp-nonOrth.dir/inline_compute-kpm-spectralOp-nonOrth.o.d -o CMakeFiles/inline_compute-kpm-spectralOp-nonOrth.dir/inline_compute-kpm-spectralOp-nonOrth.o -c /home/santiago/Documents/ICN2/Codes/linqt-3.0.0_beta/src/inline_compute-kpm-spectralOp-nonOrth.cpp

src/CMakeFiles/inline_compute-kpm-spectralOp-nonOrth.dir/inline_compute-kpm-spectralOp-nonOrth.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/inline_compute-kpm-spectralOp-nonOrth.dir/inline_compute-kpm-spectralOp-nonOrth.i"
	cd /home/santiago/Documents/ICN2/Codes/linqt-3.0.0_beta/build/src && /opt/intel/oneapi/compiler/2024.1/bin/icpx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/santiago/Documents/ICN2/Codes/linqt-3.0.0_beta/src/inline_compute-kpm-spectralOp-nonOrth.cpp > CMakeFiles/inline_compute-kpm-spectralOp-nonOrth.dir/inline_compute-kpm-spectralOp-nonOrth.i

src/CMakeFiles/inline_compute-kpm-spectralOp-nonOrth.dir/inline_compute-kpm-spectralOp-nonOrth.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/inline_compute-kpm-spectralOp-nonOrth.dir/inline_compute-kpm-spectralOp-nonOrth.s"
	cd /home/santiago/Documents/ICN2/Codes/linqt-3.0.0_beta/build/src && /opt/intel/oneapi/compiler/2024.1/bin/icpx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/santiago/Documents/ICN2/Codes/linqt-3.0.0_beta/src/inline_compute-kpm-spectralOp-nonOrth.cpp -o CMakeFiles/inline_compute-kpm-spectralOp-nonOrth.dir/inline_compute-kpm-spectralOp-nonOrth.s

# Object files for target inline_compute-kpm-spectralOp-nonOrth
inline_compute__kpm__spectralOp__nonOrth_OBJECTS = \
"CMakeFiles/inline_compute-kpm-spectralOp-nonOrth.dir/inline_compute-kpm-spectralOp-nonOrth.o"

# External object files for target inline_compute-kpm-spectralOp-nonOrth
inline_compute__kpm__spectralOp__nonOrth_EXTERNAL_OBJECTS =

inline_compute-kpm-spectralOp-nonOrth: src/CMakeFiles/inline_compute-kpm-spectralOp-nonOrth.dir/inline_compute-kpm-spectralOp-nonOrth.o
inline_compute-kpm-spectralOp-nonOrth: src/CMakeFiles/inline_compute-kpm-spectralOp-nonOrth.dir/build.make
inline_compute-kpm-spectralOp-nonOrth: lib/libkpm_lib.a
inline_compute-kpm-spectralOp-nonOrth: src/CMakeFiles/inline_compute-kpm-spectralOp-nonOrth.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/santiago/Documents/ICN2/Codes/linqt-3.0.0_beta/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../inline_compute-kpm-spectralOp-nonOrth"
	cd /home/santiago/Documents/ICN2/Codes/linqt-3.0.0_beta/build/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/inline_compute-kpm-spectralOp-nonOrth.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/inline_compute-kpm-spectralOp-nonOrth.dir/build: inline_compute-kpm-spectralOp-nonOrth
.PHONY : src/CMakeFiles/inline_compute-kpm-spectralOp-nonOrth.dir/build

src/CMakeFiles/inline_compute-kpm-spectralOp-nonOrth.dir/clean:
	cd /home/santiago/Documents/ICN2/Codes/linqt-3.0.0_beta/build/src && $(CMAKE_COMMAND) -P CMakeFiles/inline_compute-kpm-spectralOp-nonOrth.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/inline_compute-kpm-spectralOp-nonOrth.dir/clean

src/CMakeFiles/inline_compute-kpm-spectralOp-nonOrth.dir/depend:
	cd /home/santiago/Documents/ICN2/Codes/linqt-3.0.0_beta/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/santiago/Documents/ICN2/Codes/linqt-3.0.0_beta /home/santiago/Documents/ICN2/Codes/linqt-3.0.0_beta/src /home/santiago/Documents/ICN2/Codes/linqt-3.0.0_beta/build /home/santiago/Documents/ICN2/Codes/linqt-3.0.0_beta/build/src /home/santiago/Documents/ICN2/Codes/linqt-3.0.0_beta/build/src/CMakeFiles/inline_compute-kpm-spectralOp-nonOrth.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/inline_compute-kpm-spectralOp-nonOrth.dir/depend

