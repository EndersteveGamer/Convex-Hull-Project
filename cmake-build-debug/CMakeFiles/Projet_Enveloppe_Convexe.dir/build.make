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
CMAKE_SOURCE_DIR = "/mnt/c/Users/theop/OneDrive/Bureau/Université/Licence 2/Semestre 3/Algorithmique et Structure de Données/Projet Enveloppe Convexe"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/mnt/c/Users/theop/OneDrive/Bureau/Université/Licence 2/Semestre 3/Algorithmique et Structure de Données/Projet Enveloppe Convexe/cmake-build-debug"

# Include any dependencies generated for this target.
include CMakeFiles/Projet_Enveloppe_Convexe.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/Projet_Enveloppe_Convexe.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/Projet_Enveloppe_Convexe.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Projet_Enveloppe_Convexe.dir/flags.make

CMakeFiles/Projet_Enveloppe_Convexe.dir/hull.c.o: CMakeFiles/Projet_Enveloppe_Convexe.dir/flags.make
CMakeFiles/Projet_Enveloppe_Convexe.dir/hull.c.o: ../hull.c
CMakeFiles/Projet_Enveloppe_Convexe.dir/hull.c.o: CMakeFiles/Projet_Enveloppe_Convexe.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/mnt/c/Users/theop/OneDrive/Bureau/Université/Licence 2/Semestre 3/Algorithmique et Structure de Données/Projet Enveloppe Convexe/cmake-build-debug/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/Projet_Enveloppe_Convexe.dir/hull.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/Projet_Enveloppe_Convexe.dir/hull.c.o -MF CMakeFiles/Projet_Enveloppe_Convexe.dir/hull.c.o.d -o CMakeFiles/Projet_Enveloppe_Convexe.dir/hull.c.o -c "/mnt/c/Users/theop/OneDrive/Bureau/Université/Licence 2/Semestre 3/Algorithmique et Structure de Données/Projet Enveloppe Convexe/hull.c"

CMakeFiles/Projet_Enveloppe_Convexe.dir/hull.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/Projet_Enveloppe_Convexe.dir/hull.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E "/mnt/c/Users/theop/OneDrive/Bureau/Université/Licence 2/Semestre 3/Algorithmique et Structure de Données/Projet Enveloppe Convexe/hull.c" > CMakeFiles/Projet_Enveloppe_Convexe.dir/hull.c.i

CMakeFiles/Projet_Enveloppe_Convexe.dir/hull.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/Projet_Enveloppe_Convexe.dir/hull.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S "/mnt/c/Users/theop/OneDrive/Bureau/Université/Licence 2/Semestre 3/Algorithmique et Structure de Données/Projet Enveloppe Convexe/hull.c" -o CMakeFiles/Projet_Enveloppe_Convexe.dir/hull.c.s

# Object files for target Projet_Enveloppe_Convexe
Projet_Enveloppe_Convexe_OBJECTS = \
"CMakeFiles/Projet_Enveloppe_Convexe.dir/hull.c.o"

# External object files for target Projet_Enveloppe_Convexe
Projet_Enveloppe_Convexe_EXTERNAL_OBJECTS =

Projet_Enveloppe_Convexe: CMakeFiles/Projet_Enveloppe_Convexe.dir/hull.c.o
Projet_Enveloppe_Convexe: CMakeFiles/Projet_Enveloppe_Convexe.dir/build.make
Projet_Enveloppe_Convexe: CMakeFiles/Projet_Enveloppe_Convexe.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/mnt/c/Users/theop/OneDrive/Bureau/Université/Licence 2/Semestre 3/Algorithmique et Structure de Données/Projet Enveloppe Convexe/cmake-build-debug/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable Projet_Enveloppe_Convexe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Projet_Enveloppe_Convexe.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Projet_Enveloppe_Convexe.dir/build: Projet_Enveloppe_Convexe
.PHONY : CMakeFiles/Projet_Enveloppe_Convexe.dir/build

CMakeFiles/Projet_Enveloppe_Convexe.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/Projet_Enveloppe_Convexe.dir/cmake_clean.cmake
.PHONY : CMakeFiles/Projet_Enveloppe_Convexe.dir/clean

CMakeFiles/Projet_Enveloppe_Convexe.dir/depend:
	cd "/mnt/c/Users/theop/OneDrive/Bureau/Université/Licence 2/Semestre 3/Algorithmique et Structure de Données/Projet Enveloppe Convexe/cmake-build-debug" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/mnt/c/Users/theop/OneDrive/Bureau/Université/Licence 2/Semestre 3/Algorithmique et Structure de Données/Projet Enveloppe Convexe" "/mnt/c/Users/theop/OneDrive/Bureau/Université/Licence 2/Semestre 3/Algorithmique et Structure de Données/Projet Enveloppe Convexe" "/mnt/c/Users/theop/OneDrive/Bureau/Université/Licence 2/Semestre 3/Algorithmique et Structure de Données/Projet Enveloppe Convexe/cmake-build-debug" "/mnt/c/Users/theop/OneDrive/Bureau/Université/Licence 2/Semestre 3/Algorithmique et Structure de Données/Projet Enveloppe Convexe/cmake-build-debug" "/mnt/c/Users/theop/OneDrive/Bureau/Université/Licence 2/Semestre 3/Algorithmique et Structure de Données/Projet Enveloppe Convexe/cmake-build-debug/CMakeFiles/Projet_Enveloppe_Convexe.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/Projet_Enveloppe_Convexe.dir/depend

