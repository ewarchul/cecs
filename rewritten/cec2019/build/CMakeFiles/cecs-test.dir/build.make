# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_SOURCE_DIR = /home/ewarchul/cecs/rewritten/cec2019

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ewarchul/cecs/rewritten/cec2019/build

# Include any dependencies generated for this target.
include CMakeFiles/cecs-test.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/cecs-test.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/cecs-test.dir/flags.make

CMakeFiles/cecs-test.dir/test/test_basics.c.o: CMakeFiles/cecs-test.dir/flags.make
CMakeFiles/cecs-test.dir/test/test_basics.c.o: ../test/test_basics.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ewarchul/cecs/rewritten/cec2019/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/cecs-test.dir/test/test_basics.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/cecs-test.dir/test/test_basics.c.o   -c /home/ewarchul/cecs/rewritten/cec2019/test/test_basics.c

CMakeFiles/cecs-test.dir/test/test_basics.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/cecs-test.dir/test/test_basics.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/ewarchul/cecs/rewritten/cec2019/test/test_basics.c > CMakeFiles/cecs-test.dir/test/test_basics.c.i

CMakeFiles/cecs-test.dir/test/test_basics.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/cecs-test.dir/test/test_basics.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/ewarchul/cecs/rewritten/cec2019/test/test_basics.c -o CMakeFiles/cecs-test.dir/test/test_basics.c.s

CMakeFiles/cecs-test.dir/src/cec2019.c.o: CMakeFiles/cecs-test.dir/flags.make
CMakeFiles/cecs-test.dir/src/cec2019.c.o: ../src/cec2019.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ewarchul/cecs/rewritten/cec2019/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object CMakeFiles/cecs-test.dir/src/cec2019.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/cecs-test.dir/src/cec2019.c.o   -c /home/ewarchul/cecs/rewritten/cec2019/src/cec2019.c

CMakeFiles/cecs-test.dir/src/cec2019.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/cecs-test.dir/src/cec2019.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/ewarchul/cecs/rewritten/cec2019/src/cec2019.c > CMakeFiles/cecs-test.dir/src/cec2019.c.i

CMakeFiles/cecs-test.dir/src/cec2019.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/cecs-test.dir/src/cec2019.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/ewarchul/cecs/rewritten/cec2019/src/cec2019.c -o CMakeFiles/cecs-test.dir/src/cec2019.c.s

CMakeFiles/cecs-test.dir/src/cec2019_functions.c.o: CMakeFiles/cecs-test.dir/flags.make
CMakeFiles/cecs-test.dir/src/cec2019_functions.c.o: ../src/cec2019_functions.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ewarchul/cecs/rewritten/cec2019/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object CMakeFiles/cecs-test.dir/src/cec2019_functions.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/cecs-test.dir/src/cec2019_functions.c.o   -c /home/ewarchul/cecs/rewritten/cec2019/src/cec2019_functions.c

CMakeFiles/cecs-test.dir/src/cec2019_functions.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/cecs-test.dir/src/cec2019_functions.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/ewarchul/cecs/rewritten/cec2019/src/cec2019_functions.c > CMakeFiles/cecs-test.dir/src/cec2019_functions.c.i

CMakeFiles/cecs-test.dir/src/cec2019_functions.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/cecs-test.dir/src/cec2019_functions.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/ewarchul/cecs/rewritten/cec2019/src/cec2019_functions.c -o CMakeFiles/cecs-test.dir/src/cec2019_functions.c.s

CMakeFiles/cecs-test.dir/src/cec2019_interface.c.o: CMakeFiles/cecs-test.dir/flags.make
CMakeFiles/cecs-test.dir/src/cec2019_interface.c.o: ../src/cec2019_interface.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ewarchul/cecs/rewritten/cec2019/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object CMakeFiles/cecs-test.dir/src/cec2019_interface.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/cecs-test.dir/src/cec2019_interface.c.o   -c /home/ewarchul/cecs/rewritten/cec2019/src/cec2019_interface.c

CMakeFiles/cecs-test.dir/src/cec2019_interface.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/cecs-test.dir/src/cec2019_interface.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/ewarchul/cecs/rewritten/cec2019/src/cec2019_interface.c > CMakeFiles/cecs-test.dir/src/cec2019_interface.c.i

CMakeFiles/cecs-test.dir/src/cec2019_interface.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/cecs-test.dir/src/cec2019_interface.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/ewarchul/cecs/rewritten/cec2019/src/cec2019_interface.c -o CMakeFiles/cecs-test.dir/src/cec2019_interface.c.s

CMakeFiles/cecs-test.dir/src/cecs.c.o: CMakeFiles/cecs-test.dir/flags.make
CMakeFiles/cecs-test.dir/src/cecs.c.o: ../src/cecs.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ewarchul/cecs/rewritten/cec2019/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building C object CMakeFiles/cecs-test.dir/src/cecs.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/cecs-test.dir/src/cecs.c.o   -c /home/ewarchul/cecs/rewritten/cec2019/src/cecs.c

CMakeFiles/cecs-test.dir/src/cecs.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/cecs-test.dir/src/cecs.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/ewarchul/cecs/rewritten/cec2019/src/cecs.c > CMakeFiles/cecs-test.dir/src/cecs.c.i

CMakeFiles/cecs-test.dir/src/cecs.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/cecs-test.dir/src/cecs.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/ewarchul/cecs/rewritten/cec2019/src/cecs.c -o CMakeFiles/cecs-test.dir/src/cecs.c.s

# Object files for target cecs-test
cecs__test_OBJECTS = \
"CMakeFiles/cecs-test.dir/test/test_basics.c.o" \
"CMakeFiles/cecs-test.dir/src/cec2019.c.o" \
"CMakeFiles/cecs-test.dir/src/cec2019_functions.c.o" \
"CMakeFiles/cecs-test.dir/src/cec2019_interface.c.o" \
"CMakeFiles/cecs-test.dir/src/cecs.c.o"

# External object files for target cecs-test
cecs__test_EXTERNAL_OBJECTS =

../bin/cecs-test: CMakeFiles/cecs-test.dir/test/test_basics.c.o
../bin/cecs-test: CMakeFiles/cecs-test.dir/src/cec2019.c.o
../bin/cecs-test: CMakeFiles/cecs-test.dir/src/cec2019_functions.c.o
../bin/cecs-test: CMakeFiles/cecs-test.dir/src/cec2019_interface.c.o
../bin/cecs-test: CMakeFiles/cecs-test.dir/src/cecs.c.o
../bin/cecs-test: CMakeFiles/cecs-test.dir/build.make
../bin/cecs-test: /usr/local/lib/libcmocka.so
../bin/cecs-test: CMakeFiles/cecs-test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/ewarchul/cecs/rewritten/cec2019/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking C executable ../bin/cecs-test"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/cecs-test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/cecs-test.dir/build: ../bin/cecs-test

.PHONY : CMakeFiles/cecs-test.dir/build

CMakeFiles/cecs-test.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/cecs-test.dir/cmake_clean.cmake
.PHONY : CMakeFiles/cecs-test.dir/clean

CMakeFiles/cecs-test.dir/depend:
	cd /home/ewarchul/cecs/rewritten/cec2019/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ewarchul/cecs/rewritten/cec2019 /home/ewarchul/cecs/rewritten/cec2019 /home/ewarchul/cecs/rewritten/cec2019/build /home/ewarchul/cecs/rewritten/cec2019/build /home/ewarchul/cecs/rewritten/cec2019/build/CMakeFiles/cecs-test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/cecs-test.dir/depend

