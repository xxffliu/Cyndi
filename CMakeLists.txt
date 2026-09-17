cmake_minimum_required(VERSION 3.15)
project(Cyndi VERSION 4.0 LANGUAGES CXX)

# Cyndi: multi-objective evolutionary generation of bioactive conformers.
# Liu et al., BMC Bioinformatics 10, 101 (2009).
#
#   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
#   cmake --build build -j
#   ctest --test-dir build --output-on-failure
#
# On Windows with the MinGW-w64 toolchain in tools/mingw64, use the preset:
#
#   cmake --preset mingw-release
#   cmake --build --preset mingw-release
#   ctest --preset mingw-release
#
# Options:
#   -DCYNDI_BUILD_TESTS=ON       build and register the smoke tests (default ON)
#   -DCYNDI_WERROR=ON            turn all warnings into errors (default OFF -- see below)
#   -DCYNDI_WERROR_REORDER=ON    turn out-of-order member initialisers into errors
#                                (default OFF -- 167 of them remain, see below)
#   -DCYNDI_TORSION_LINEAR=ON    build the pre-v1 linear torsion encoding (default OFF)
#   -DCYNDI_SANITIZE=address|thread|undefined
#                                build with the given sanitizer (default: none)
#   -DCYNDI_STATIC_RUNTIME=ON    link libgcc/libstdc++ statically (default ON on MinGW)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
	set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(CYNDI_BUILD_TESTS     "Build the smoke and regression tests" ON)
option(CYNDI_WERROR          "Treat warnings as errors" OFF)
option(CYNDI_WERROR_REORDER  "Treat out-of-order member initialisers as errors" OFF)
option(CYNDI_TORSION_LINEAR  "Use the original linear (non-periodic) torsion encoding" OFF)
if(MINGW)
	option(CYNDI_STATIC_RUNTIME "Link the GCC runtime statically" ON)
else()
	option(CYNDI_STATIC_RUNTIME "Link the GCC runtime statically" OFF)
endif()
set(CYNDI_SANITIZE "" CACHE STRING "Sanitizer to build with: address, thread, undefined or empty")

# ---------------------------------------------------------------------------
# warnings
#
# CYNDI_WERROR defaults to OFF, unlike SHAFTS. This tree still emits ~750
# -Wall -Wextra warnings; they have not been worked through, so -Werror would
# simply fail the build. The two categories worth attention when someone does
# get to them are -Wmisleading-indentation and -Wdangling-else, concentrated in
# the MMFF94 parameter-assignment rule cascades (MMFF94Torsion.cpp around lines
# 270-355 in particular). They were read in 2026-09 and are cosmetic -- the
# `if (!found_rule)` guards are sequential, just progressively indented -- but
# that is exactly the shape a real bug hides in, so do not let the count grow.
#
# -Wreorder is worth its own switch. Member initialisers that do not follow
# declaration order are how an uninitialised member hides -- a constructor that
# LOOKS like it sets everything, in an order the compiler silently ignores.
# SHAFTS turns it into an error unconditionally; this tree cannot yet, because
# it has 167 of them across 18 constructors, 108 of those in the single
# MOGAParam constructor in include/ParamInput.h. They are all harmless today
# (every initialiser is a literal, so declaration order gives the same result),
# but that is luck, not design. Cleaning them up and flipping
# CYNDI_WERROR_REORDER to ON is a good small task for whoever is next in here.
# ---------------------------------------------------------------------------
add_library(cyndi_warnings INTERFACE)
if(MSVC)
	target_compile_options(cyndi_warnings INTERFACE /W3)
	if(CYNDI_WERROR)
		target_compile_options(cyndi_warnings INTERFACE /WX)
	endif()
else()
	target_compile_options(cyndi_warnings INTERFACE -Wall -Wextra -Wreorder)
	if(CYNDI_WERROR_REORDER)
		target_compile_options(cyndi_warnings INTERFACE -Werror=reorder)
	endif()
	if(CYNDI_WERROR)
		target_compile_options(cyndi_warnings INTERFACE -Werror)
	endif()
endif()

add_library(cyndi_sanitizers INTERFACE)
if(CYNDI_SANITIZE)
	if(MSVC)
		message(WARNING "CYNDI_SANITIZE is ignored with MSVC")
	else()
		target_compile_options(cyndi_sanitizers INTERFACE -fsanitize=${CYNDI_SANITIZE} -fno-omit-frame-pointer -g)
		target_link_options(cyndi_sanitizers INTERFACE -fsanitize=${CYNDI_SANITIZE})
	endif()
endif()

# ---------------------------------------------------------------------------
# MinGW sysroot headers
#
# The MinGW-w64 build in tools/mingw64 does not put its own <prefix>/include on
# the default search path (this is what build.sh worked around with an explicit
# -idirafter). Detect it from the compiler location and add it the same way:
# -idirafter, so it is searched last and never shadows this project's headers.
# Harmless on a correctly configured toolchain, since the directory is then
# already on the path.
# ---------------------------------------------------------------------------
add_library(cyndi_platform INTERFACE)
if(MINGW)
	get_filename_component(_cyndi_mingw_bin "${CMAKE_CXX_COMPILER}" DIRECTORY)
	get_filename_component(_cyndi_mingw_root "${_cyndi_mingw_bin}" DIRECTORY)
	if(EXISTS "${_cyndi_mingw_root}/include")
		target_compile_options(cyndi_platform INTERFACE "-idirafter" "${_cyndi_mingw_root}/include")
		message(STATUS "MinGW sysroot headers: ${_cyndi_mingw_root}/include (added with -idirafter)")
	endif()
endif()
if(CYNDI_STATIC_RUNTIME AND NOT MSVC)
	# Without this the executable needs libstdc++-6.dll, libgcc_s_seh-1.dll and
	# libwinpthread-1.dll next to it -- which is why build.sh copied them into
	# the output directory. Linking them in removes the whole problem.
	target_link_options(cyndi_platform INTERFACE -static-libgcc -static-libstdc++)
	if(MINGW)
		target_link_options(cyndi_platform INTERFACE -static)
	endif()
endif()

# ---------------------------------------------------------------------------
# core library: every source except the one holding main()
# ---------------------------------------------------------------------------
set(CYNDI_CORE_SOURCES
	src/Atom.cpp
	src/Bond.cpp
	src/CGMinimizer.cpp
	src/Cluster.cpp
	src/ConGen.cpp
	src/ConicLineSearch.cpp
	src/EnergyMinimizer.cpp
	src/FFComponent.cpp
	src/FFParameter.cpp
	src/ForceField.cpp
	src/Fragment.cpp
	src/Gradient.cpp
	src/MMFF94.cpp
	src/MMFF94Bend.cpp
	src/MMFF94Electrostatic.cpp
	src/MMFF94OOP.cpp
	src/MMFF94Stretch.cpp
	src/MMFF94StretchBend.cpp
	src/MMFF94Torsion.cpp
	src/MMFF94VDW.cpp
	src/MOL2IO.cpp
	src/Matrix3x3.cpp
	src/Mol.cpp
	src/ParamInput.cpp
	src/SDMinimizer.cpp
	src/TAFF.cpp
	src/TAFFBend.cpp
	src/TAFFElectrostatic.cpp
	src/TAFFOOP.cpp
	src/TAFFStretch.cpp
	src/TAFFTorsion.cpp
	src/TAFFVDW.cpp
	src/Vector3.cpp
	src/random.cpp
)

add_library(cyndi_core STATIC ${CYNDI_CORE_SOURCES})
target_include_directories(cyndi_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(cyndi_core PUBLIC cyndi_warnings cyndi_sanitizers cyndi_platform)

# TORSION_LINEAR switches crossover and mutation back to clamping at +/-180 deg
# instead of wrapping around the circle. It changes results, so it is a whole-
# project property rather than a per-file one.
if(CYNDI_TORSION_LINEAR)
	target_compile_definitions(cyndi_core PUBLIC TORSION_LINEAR)
	message(STATUS "Torsion encoding: LINEAR (pre-v1 behaviour)")
else()
	message(STATUS "Torsion encoding: circular (default)")
endif()

# ---------------------------------------------------------------------------
# executable
# ---------------------------------------------------------------------------
add_executable(cyndi src/Cyndi.cpp)
target_link_libraries(cyndi PRIVATE cyndi_core)

# MMFF94.parm and TAFF.parm are opened by bare filename (see
# DEFAULT_MMFF94_PARAM_FILE in include/MMFF94.h), i.e. relative to the working
# directory, not to the executable. Put copies next to the binary so that
# running it from the build directory just works; a real fix would be an
# explicit -ffparam option.
set(CYNDI_RUNTIME_DIR $<TARGET_FILE_DIR:cyndi>)
set(CYNDI_PARM_FILES "")
foreach(_parm MMFF94.parm TAFF.parm)
	if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/Cyndi/${_parm}")
		list(APPEND CYNDI_PARM_FILES "${CMAKE_CURRENT_SOURCE_DIR}/Cyndi/${_parm}")
		add_custom_command(TARGET cyndi POST_BUILD
			COMMAND ${CMAKE_COMMAND} -E copy_if_different
				"${CMAKE_CURRENT_SOURCE_DIR}/Cyndi/${_parm}" "${CYNDI_RUNTIME_DIR}/${_parm}"
			COMMENT "Staging ${_parm} next to the executable")
	else()
		message(WARNING "Cyndi/${_parm} is missing -- the executable will not run "
		                "until it is on the working directory")
	endif()
endforeach()

# ---------------------------------------------------------------------------
# tests
#
# There is no unit-test suite in this tree. What is registered here are the
# end-to-end invariants that the 2026-09 review actually established, so that
# they cannot regress silently:
#
#   runs_single_molecule   the program completes on a small ligand
#   chunk_reproducibility  a chunked run (-startidx/-maxmols) reproduces the
#                          same molecule from a whole-file run byte for byte.
#                          This is KNOWN_ISSUES #1 and it has broken twice --
#                          once in GlobalIdx, once in MOL2IO stream state.
#   seed_wraparound        molecules past index 580 do not segfault
#                          (KNOWN_ISSUES #4)
#   rmsd_evaluator         bench/rmsd_core.py agrees with a reference rotation
#   bad_option_rejected    an unknown CLI option exits non-zero
# ---------------------------------------------------------------------------
if(CYNDI_BUILD_TESTS)
	enable_testing()

	set(CYNDI_TEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/tests-run)
	file(MAKE_DIRECTORY ${CYNDI_TEST_DIR})
	foreach(_parm ${CYNDI_PARM_FILES})
		configure_file(${_parm} ${CYNDI_TEST_DIR} COPYONLY)
	endforeach()

	# A deliberately tiny search so the suite stays a few seconds, not minutes.
	file(WRITE ${CYNDI_TEST_DIR}/test_param.in
"MOGA_Max_Conformers\t\t20
MOGA_Num_Objectives\t\t3
MOGA_Max_Generation\t\t20
MOGA_Population_Size\t\t30
MOGA_Max_Run\t\t1
MOGA_Crossover_Probability\t0.85
MOGA_Mutation_Probability\t0.4
MOGA_SBX\t\t\t15 20
MOGA_Epsilon_Quaternion\t\t2 0.2 0.1 2
MOGA_Keep_Input_Conformer[Y/N]\tN
MOGA_Optimize_Conformer[Y/N]\tY
MOGA_Energy_Cutoff\t\t60
MOGA_RMSD_Scale_Factor\t\t0.2
MOGA_Max_Opt_Iteration\t\t20
MOGA_Max_Opt_Gradient\t\t0.5
MOGA_Random_Seed\t\t0.42
")

	add_test(NAME runs_single_molecule
	         COMMAND cyndi -input ${CMAKE_CURRENT_SOURCE_DIR}/Cyndi/3ert.mol2
	                       -output ${CYNDI_TEST_DIR}/3ert_out.mol2
	                       -parm ${CYNDI_TEST_DIR}/test_param.in)
	set_tests_properties(runs_single_molecule PROPERTIES
		WORKING_DIRECTORY ${CYNDI_TEST_DIR}
		PASS_REGULAR_EXPRESSION "Number of Molecules Failed: 0")

	add_test(NAME bad_option_rejected COMMAND cyndi -input x.mol2 -nonsense)
	set_tests_properties(bad_option_rejected PROPERTIES
		WORKING_DIRECTORY ${CYNDI_TEST_DIR} WILL_FAIL TRUE)

	find_package(Python3 COMPONENTS Interpreter)
	if(Python3_Interpreter_FOUND)
		# These three need to drive the binary and compare files, so they are
		# python rather than raw ctest commands. Each skips itself (exit 77,
		# SKIP_RETURN_CODE below) when its input data is not in the tree.
		add_test(NAME chunk_reproducibility
		         COMMAND ${Python3_EXECUTABLE}
		                 ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_chunk_reproducibility.py
		                 $<TARGET_FILE:cyndi>
		                 ${CMAKE_CURRENT_SOURCE_DIR}/bench/329_test_set.mol2
		                 ${CYNDI_TEST_DIR}/test_param.in)
		set_tests_properties(chunk_reproducibility PROPERTIES
			WORKING_DIRECTORY ${CYNDI_TEST_DIR} SKIP_RETURN_CODE 77 TIMEOUT 1800)

		add_test(NAME seed_wraparound
		         COMMAND ${Python3_EXECUTABLE}
		                 ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_seed_wraparound.py
		                 $<TARGET_FILE:cyndi>
		                 ${CMAKE_CURRENT_SOURCE_DIR}/bench/329_test_set.mol2
		                 ${CYNDI_TEST_DIR}/test_param.in)
		set_tests_properties(seed_wraparound PROPERTIES
			WORKING_DIRECTORY ${CYNDI_TEST_DIR} SKIP_RETURN_CODE 77 TIMEOUT 900)

		add_test(NAME rmsd_evaluator
		         COMMAND ${Python3_EXECUTABLE}
		                 ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_rmsd_core.py)
		set_tests_properties(rmsd_evaluator PROPERTIES SKIP_RETURN_CODE 77)
	else()
		message(STATUS "Python3 not found -- chunk_reproducibility, seed_wraparound "
		               "and rmsd_evaluator will not be registered")
	endif()
endif()

install(TARGETS cyndi RUNTIME DESTINATION bin)
if(CYNDI_PARM_FILES)
	install(FILES ${CYNDI_PARM_FILES} DESTINATION bin)
endif()
install(FILES parameter/CyndiParam.in DESTINATION share/cyndi)
