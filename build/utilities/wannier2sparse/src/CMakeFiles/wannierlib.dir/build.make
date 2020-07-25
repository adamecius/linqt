# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
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
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /data/jgarcia/codes/linqt

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /data/jgarcia/codes/linqt/build

# Include any dependencies generated for this target.
include utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/depend.make

# Include the progress variables for this target.
include utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/progress.make

# Include the compile flags for this target's objects.
include utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/flags.make

utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/wannier_parser.cpp.o: utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/flags.make
utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/wannier_parser.cpp.o: ../utilities/wannier2sparse/src/wannier_parser.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/data/jgarcia/codes/linqt/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/wannier_parser.cpp.o"
	cd /data/jgarcia/codes/linqt/build/utilities/wannier2sparse/src && /home/ICN2/jgarcia/.local/intel/compilers_and_libraries_2019.1.144/linux/bin/intel64/icpc   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/wannierlib.dir/wannier_parser.cpp.o -c /data/jgarcia/codes/linqt/utilities/wannier2sparse/src/wannier_parser.cpp

utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/wannier_parser.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/wannierlib.dir/wannier_parser.cpp.i"
	cd /data/jgarcia/codes/linqt/build/utilities/wannier2sparse/src && /home/ICN2/jgarcia/.local/intel/compilers_and_libraries_2019.1.144/linux/bin/intel64/icpc  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /data/jgarcia/codes/linqt/utilities/wannier2sparse/src/wannier_parser.cpp > CMakeFiles/wannierlib.dir/wannier_parser.cpp.i

utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/wannier_parser.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/wannierlib.dir/wannier_parser.cpp.s"
	cd /data/jgarcia/codes/linqt/build/utilities/wannier2sparse/src && /home/ICN2/jgarcia/.local/intel/compilers_and_libraries_2019.1.144/linux/bin/intel64/icpc  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /data/jgarcia/codes/linqt/utilities/wannier2sparse/src/wannier_parser.cpp -o CMakeFiles/wannierlib.dir/wannier_parser.cpp.s

utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/wannier_parser.cpp.o.requires:

.PHONY : utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/wannier_parser.cpp.o.requires

utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/wannier_parser.cpp.o.provides: utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/wannier_parser.cpp.o.requires
	$(MAKE) -f utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/build.make utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/wannier_parser.cpp.o.provides.build
.PHONY : utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/wannier_parser.cpp.o.provides

utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/wannier_parser.cpp.o.provides.build: utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/wannier_parser.cpp.o


utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/hopping_list.cpp.o: utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/flags.make
utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/hopping_list.cpp.o: ../utilities/wannier2sparse/src/hopping_list.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/data/jgarcia/codes/linqt/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/hopping_list.cpp.o"
	cd /data/jgarcia/codes/linqt/build/utilities/wannier2sparse/src && /home/ICN2/jgarcia/.local/intel/compilers_and_libraries_2019.1.144/linux/bin/intel64/icpc   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/wannierlib.dir/hopping_list.cpp.o -c /data/jgarcia/codes/linqt/utilities/wannier2sparse/src/hopping_list.cpp

utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/hopping_list.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/wannierlib.dir/hopping_list.cpp.i"
	cd /data/jgarcia/codes/linqt/build/utilities/wannier2sparse/src && /home/ICN2/jgarcia/.local/intel/compilers_and_libraries_2019.1.144/linux/bin/intel64/icpc  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /data/jgarcia/codes/linqt/utilities/wannier2sparse/src/hopping_list.cpp > CMakeFiles/wannierlib.dir/hopping_list.cpp.i

utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/hopping_list.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/wannierlib.dir/hopping_list.cpp.s"
	cd /data/jgarcia/codes/linqt/build/utilities/wannier2sparse/src && /home/ICN2/jgarcia/.local/intel/compilers_and_libraries_2019.1.144/linux/bin/intel64/icpc  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /data/jgarcia/codes/linqt/utilities/wannier2sparse/src/hopping_list.cpp -o CMakeFiles/wannierlib.dir/hopping_list.cpp.s

utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/hopping_list.cpp.o.requires:

.PHONY : utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/hopping_list.cpp.o.requires

utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/hopping_list.cpp.o.provides: utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/hopping_list.cpp.o.requires
	$(MAKE) -f utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/build.make utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/hopping_list.cpp.o.provides.build
.PHONY : utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/hopping_list.cpp.o.provides

utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/hopping_list.cpp.o.provides.build: utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/hopping_list.cpp.o


utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/tbmodel.cpp.o: utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/flags.make
utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/tbmodel.cpp.o: ../utilities/wannier2sparse/src/tbmodel.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/data/jgarcia/codes/linqt/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/tbmodel.cpp.o"
	cd /data/jgarcia/codes/linqt/build/utilities/wannier2sparse/src && /home/ICN2/jgarcia/.local/intel/compilers_and_libraries_2019.1.144/linux/bin/intel64/icpc   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/wannierlib.dir/tbmodel.cpp.o -c /data/jgarcia/codes/linqt/utilities/wannier2sparse/src/tbmodel.cpp

utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/tbmodel.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/wannierlib.dir/tbmodel.cpp.i"
	cd /data/jgarcia/codes/linqt/build/utilities/wannier2sparse/src && /home/ICN2/jgarcia/.local/intel/compilers_and_libraries_2019.1.144/linux/bin/intel64/icpc  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /data/jgarcia/codes/linqt/utilities/wannier2sparse/src/tbmodel.cpp > CMakeFiles/wannierlib.dir/tbmodel.cpp.i

utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/tbmodel.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/wannierlib.dir/tbmodel.cpp.s"
	cd /data/jgarcia/codes/linqt/build/utilities/wannier2sparse/src && /home/ICN2/jgarcia/.local/intel/compilers_and_libraries_2019.1.144/linux/bin/intel64/icpc  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /data/jgarcia/codes/linqt/utilities/wannier2sparse/src/tbmodel.cpp -o CMakeFiles/wannierlib.dir/tbmodel.cpp.s

utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/tbmodel.cpp.o.requires:

.PHONY : utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/tbmodel.cpp.o.requires

utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/tbmodel.cpp.o.provides: utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/tbmodel.cpp.o.requires
	$(MAKE) -f utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/build.make utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/tbmodel.cpp.o.provides.build
.PHONY : utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/tbmodel.cpp.o.provides

utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/tbmodel.cpp.o.provides.build: utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/tbmodel.cpp.o


# Object files for target wannierlib
wannierlib_OBJECTS = \
"CMakeFiles/wannierlib.dir/wannier_parser.cpp.o" \
"CMakeFiles/wannierlib.dir/hopping_list.cpp.o" \
"CMakeFiles/wannierlib.dir/tbmodel.cpp.o"

# External object files for target wannierlib
wannierlib_EXTERNAL_OBJECTS =

lib/libwannierlib.a: utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/wannier_parser.cpp.o
lib/libwannierlib.a: utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/hopping_list.cpp.o
lib/libwannierlib.a: utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/tbmodel.cpp.o
lib/libwannierlib.a: utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/build.make
lib/libwannierlib.a: utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/data/jgarcia/codes/linqt/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX static library ../../../lib/libwannierlib.a"
	cd /data/jgarcia/codes/linqt/build/utilities/wannier2sparse/src && $(CMAKE_COMMAND) -P CMakeFiles/wannierlib.dir/cmake_clean_target.cmake
	cd /data/jgarcia/codes/linqt/build/utilities/wannier2sparse/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/wannierlib.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/build: lib/libwannierlib.a

.PHONY : utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/build

utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/requires: utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/wannier_parser.cpp.o.requires
utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/requires: utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/hopping_list.cpp.o.requires
utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/requires: utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/tbmodel.cpp.o.requires

.PHONY : utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/requires

utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/clean:
	cd /data/jgarcia/codes/linqt/build/utilities/wannier2sparse/src && $(CMAKE_COMMAND) -P CMakeFiles/wannierlib.dir/cmake_clean.cmake
.PHONY : utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/clean

utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/depend:
	cd /data/jgarcia/codes/linqt/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /data/jgarcia/codes/linqt /data/jgarcia/codes/linqt/utilities/wannier2sparse/src /data/jgarcia/codes/linqt/build /data/jgarcia/codes/linqt/build/utilities/wannier2sparse/src /data/jgarcia/codes/linqt/build/utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : utilities/wannier2sparse/src/CMakeFiles/wannierlib.dir/depend
