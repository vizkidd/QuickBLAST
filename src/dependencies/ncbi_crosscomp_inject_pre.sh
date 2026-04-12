#!/bin/bash

#Get script directory from mingw shell for sysroot
SCRIPTPATH="$(
    cd -- "$(dirname "$0")" >/dev/null 2>&1
    pwd -P
)"

FOLDER_MODE="-m" #for cygpath
NCBI_DIR=$(cygpath $FOLDER_MODE $(pwd))
SCRIPTPATH=$(cygpath $FOLDER_MODE $SCRIPTPATH)
DEP_BUILD_DIR="BUILD"
RTOOLS_DIR=$(cygpath $FOLDER_MODE $RTOOLS_DIR)

OPENMP_LIB_FLAG=""
if [[ ! -z $SHLIB_OPENMP_CXXFLAGS ]]; then
    OPENMP_LIB_FLAG="-fopenmp -lomp " #-lgomp.dll -lgomp"
    # SHLIB_OPENMP_CXXFLAGS="-D_OPENMP $SHLIB_OPENMP_CXXFLAGS"
else
    echo "ERROR: QuickBLAST and ncbi-c++ toolkit requires OpenMP support. Exiting"
    exit 255
fi

echo "PRE CONFIGURE INJECT..."
echo "BUILD DIR: $NCBI_DIR"
echo "SCRIPTPATH: $SCRIPTPATH"
echo "RTOOLS DIR: $RTOOLS_DIR"
echo "OpenMP Flags : $SHLIB_OPENMP_CXXFLAGS"

# echo "Creating GCC Scrubber wrapper..."
# mkdir -p fake_bin

# # 1. Create the bash scrubbers (which handle the logic safely)
# cat << 'EOF' > fake_bin/gcc_scrubber.sh
# #!/bin/bash
# REAL_GCC="C:/rtools45/x86_64-w64-mingw32.static.posix/bin/gcc.exe"
# if [[ "$*" == *"-print-"* || "$*" == *"-v"* ]]; then
#     t_out=$(mktemp)
#     t_err=$(mktemp)
#     "$REAL_GCC" "$@" > "$t_out" 2> "$t_err"
#     ret=$?
#     tr '\\' '/' < "$t_out"
#     tr '\\' '/' < "$t_err" >&2
#     rm -f "$t_out" "$t_err"
#     exit $ret
# else
#     exec "$REAL_GCC" "$@"
# fi
# EOF

# cat << 'EOF' > fake_bin/g++_scrubber.sh
# #!/bin/bash
# REAL_GCC="C:/rtools45/x86_64-w64-mingw32.static.posix/bin/g++.exe"
# if [[ "$*" == *"-print-"* || "$*" == *"-v"* ]]; then
#     t_out=$(mktemp)
#     t_err=$(mktemp)
#     "$REAL_GCC" "$@" > "$t_out" 2> "$t_err"
#     ret=$?
#     tr '\\' '/' < "$t_out"
#     tr '\\' '/' < "$t_err" >&2
#     rm -f "$t_out" "$t_err"
#     exit $ret
# else
#     exec "$REAL_GCC" "$@"
# fi
# EOF

# chmod +x fake_bin/gcc_scrubber.sh
# chmod +x fake_bin/g++_scrubber.sh

# # 2. Create the .bat files that Windows/CMake can actually execute
# cat << 'EOF' > fake_bin/gcc.bat
# @echo off
# bash "%~dp0gcc_scrubber.sh" %*
# EOF

# cat << 'EOF' > fake_bin/g++.bat
# @echo off
# bash "%~dp0g++_scrubber.sh" %*
# EOF

cd $NCBI_DIR

# # ==========================================
# # QUICKBLAST INJECTION: NEUTRALIZE CYGWIN
# # Append absolute overrides to the Toolkit's OS check script.
# # This forces the Toolkit to evaluate the build as Native Windows 
# # BEFORE it processes project requirements.
# # ==========================================
# echo "Injecting Native Windows OS Overrides into CMakeChecks..."
# cat << 'EOF' >> src/build-system/cmake/CMakeChecks.cmake

# # --- QUICKBLAST FORCE OVERRIDES ---
# set(NCBI_OS_CYGWIN OFF)
# set(NCBI_OS_MSYS OFF)
# set(NCBI_OS_MSWIN ON)
# set(CYGWIN OFF)
# set(MSYS OFF)
# set(WIN32 ON)
# set(MINGW ON)
# set(HOST_OS "windows")
# set(HOST "-mingw32")
# # ----------------------------------
# EOF

#export CC="/usr/bin/gcc"
#export CXX="/usr/bin/g++"

# #WE can ASSUME C:/rtools45/x86_64-w64-mingw32.static.posix/bin/ IS IN PATH
# # export CC="$RTOOLS_DIR/x86_64-w64-mingw32.static.posix/bin/gcc.exe"
# # export CXX="$RTOOLS_DIR/x86_64-w64-mingw32.static.posix/bin/g++.exe"
# # export C_LD="$RTOOLS_DIR/x86_64-w64-mingw32.static.posix/bin/ld.exe"
# # export CXX_LD="$RTOOLS_DIR/x86_64-w64-mingw32.static.posix/bin/ld.exe"
# # export CC="gcc.bat"
# # export CXX="g++.bat"
# export CC="gcc.exe"
# export CXX="g++.exe"
# # export C_LD="ld.exe"
# # export CXX_LD="ld.exe"

# # Generate guaranteed forward-slashed paths using cygpath -m
RTOOLS_BIN=$(cygpath -m "$RTOOLS_DIR/x86_64-w64-mingw32.static.posix/bin")
export CC="$RTOOLS_BIN/gcc.exe"
export CXX="$RTOOLS_BIN/g++.exe"
export C_LD="$RTOOLS_BIN/ld.exe"
export CXX_LD="$RTOOLS_BIN/ld.exe"
export AR_PATH="$RTOOLS_BIN/ar.exe"
export NM_PATH="$RTOOLS_BIN/nm.exe"
export RANLIB_PATH="$RTOOLS_BIN/ranlib.exe"

COMPILER_VER=$($CC -dumpversion | tr -d '[:space:]')
COMPILER_VER_NOPUNCT=$($CC -dumpversion | tr -d '[:space:]' | tr -d '[:punct:]')

toolchain_file=$(src/build-system/cmake/toolchains/cmkTool.sh GCC $COMPILER_VER | grep "in" | awk -F' ' '{print $NF}')

# # =========================================================================
# # THE GCC SCRUBBER (Neutralizes the CMake \x ABI Bug)
# # =========================================================================
# echo "Creating GCC Scrubber wrappers..."
# mkdir -p fake_bin
# export FAKE_BIN_DIR="$(pwd)/fake_bin"

# # Create a fake gcc that translates \ to / during CMake ABI checks
# cat << 'EOF' > fake_bin/gcc
# #!/bin/bash
# REAL_GCC="C:/rtools45/x86_64-w64-mingw32.static.posix/bin/gcc.exe"
# if [[ "$*" == *"-print-"* || "$*" == *"-v "* || "$*" == *"-dump"* ]]; then
#     "$REAL_GCC" "$@" 2>&1 | tr '\\' '/'
# else
#     exec "$REAL_GCC" "$@"
# fi
# EOF

# # Create a fake g++ that translates \ to / during CMake ABI checks
# cat << 'EOF' > fake_bin/g++
# #!/bin/bash
# REAL_GXX="C:/rtools45/x86_64-w64-mingw32.static.posix/bin/g++.exe"
# if [[ "$*" == *"-print-"* || "$*" == *"-v "* || "$*" == *"-dump"* ]]; then
#     "$REAL_GXX" "$@" 2>&1 | tr '\\' '/'
# else
#     exec "$REAL_GXX" "$@"
# fi
# EOF

# chmod +x fake_bin/gcc
# chmod +x fake_bin/g++

# export CC="$FAKE_BIN_DIR/gcc"
# export CXX="$FAKE_BIN_DIR/g++"
# export C_LD="$RTOOLS_BIN/ld.exe"
# export CXX_LD="$RTOOLS_BIN/ld.exe"

# export CC="gcc"
# export CXX="g++"
# export C_LD="ld"
# export CXX_LD="ld"

# # Generate the toolchain file as usual
# COMPILER_VER=$("$CC" -dumpversion | tr -d '[:space:]')
# COMPILER_VER_NOPUNCT=$("$CC" -dumpversion | tr -d '[:space:]' | tr -d '[:punct:]')
# toolchain_file=$(src/build-system/cmake/toolchains/cmkTool.sh GCC $COMPILER_VER | grep "in" | awk -F' ' '{print $NF}')

if [[ -z $toolchain_file ]]; then
    echo "EXITING. COULD NOT GUESS TOOLCHAIN"
    exit 0
fi

SHELL_TYPE=$(echo $MSYSTEM_CARCH | tr "[:upper:]" "[:lower:]")
echo "TOOLCHAIN FILE : $toolchain_file"
echo "CC : $CC"
echo "CXX : $CXX"
echo "COMPILER VER : $COMPILER_VER"

#ln -sf $CC $(realpath $CC)-$COMPILER_VER
#ln -sf $CXX $(realpath $CXX)-$COMPILER_VER

# echo "Bypassing OS gates: Locating the hidden CMake config template..."

# UNIX_TEMPLATE=$(find src/build-system -type f \( -name "*conf*.h.in" -o -name "config*.h.in" \) | head -n 1)

# if [ -z "$UNIX_TEMPLATE" ]; then
#     echo "FATAL: Could not find ANY configuration template (.h.in) in src/build-system!"
#     exit 1
# fi

# echo "Success! Found the actual template at: $UNIX_TEMPLATE"

# # 1. FIX: Path error (/DNDEBUG to -DNDEBUG)
# echo "Patching MSVC flags for MinGW GCC compatibility..."
# find src/build-system -type f -name "*.cmake" -exec sed -i 's|/DNDEBUG|-DNDEBUG|g' {} +

# # 2. FIX: ssize_t clash
# echo "Patching ssize_t definitions in the template..."
# sed -i 's/typedef __int64 ssize_t;/\/\/ ssize_t already defined by MinGW/g' "$UNIX_TEMPLATE"

# # 3. FIX: Stop CMake from defining NCBI_COMPILER_MSVC for MinGW
# echo "Stripping false MSVC compiler definitions..."
# sed -i 's/#define  NCBI_COMPILER_MSVC 1/\/* #undef NCBI_COMPILER_MSVC *\//g' "$UNIX_TEMPLATE"
# sed -i 's/#cmakedefine NCBI_COMPILER_MSVC 1/\/* #undef NCBI_COMPILER_MSVC *\//g' "$UNIX_TEMPLATE"

# # 4. FIX: Patch ncbistre.hpp so it stops trying to use Microsoft _Openprot extensions
# echo "Patching ncbistre.hpp to use standard GNU libstdc++ file streams..."
# sed -i 's/#if defined(NCBI_OS_MSWIN)/#if defined(NCBI_OS_MSWIN) \&\& defined(_MSC_VER)/g' include/corelib/ncbistre.hpp
# sed -i 's/#ifdef NCBI_OS_MSWIN/#if defined(NCBI_OS_MSWIN) \&\& defined(_MSC_VER)/g' include/corelib/ncbistre.hpp

# # 5. FIX: Destroy the hardcoded MSVC requirement gate
# echo "Removing hardcoded MSVC restriction..."
# sed -i 's/#    error "This toolkit is not buildable with a compiler other than MSVC."/\/\/ Bypassed MSVC requirement for MinGW/g' include/common/ncbi_export.h

# # 6. FIX: Resolve the constructor type collision in CProcess
# echo "Patching CProcess constructor collision..."
# sed -i 's/CProcess(TProcessHandle process, EType type = eHandle);/\/\/ CProcess(TProcessHandle process, EType type = eHandle);/g' include/corelib/ncbi_process.hpp

# # # 7. FIX: Cure the OS identity crisis
# # echo "Stripping false Unix/Linux OS macros..."
# # sed -i 's/#define NCBI_OS_UNIX 1/\/* #undef NCBI_OS_UNIX *\//g' "$UNIX_TEMPLATE"
# # sed -i 's/#cmakedefine NCBI_OS_UNIX 1/\/* #undef NCBI_OS_UNIX *\//g' "$UNIX_TEMPLATE"
# # sed -i 's/#define NCBI_OS_LINUX 1/\/* #undef NCBI_OS_LINUX *\//g' "$UNIX_TEMPLATE"
# # sed -i 's/#cmakedefine NCBI_OS_LINUX 1/\/* #undef NCBI_OS_LINUX *\//g' "$UNIX_TEMPLATE"
# # sed -i 's/#define NCBI_OS "Linux"/#define NCBI_OS "MSWIN"/g' "$UNIX_TEMPLATE"

# # # Safely patch sys/resource.h in ncbi_process.cpp just in case
# # sed -i 's|#  include <sys/resource.h>|#if defined(HAVE_SYS_RESOURCE_H)\n#  include <sys/resource.h>\n#endif|g' src/corelib/ncbi_process.cpp
# # sed -i 's|#  include <sys/time.h>|#if defined(HAVE_SYS_TIME_H)\n#  include <sys/time.h>\n#endif|g' src/corelib/ncbi_process.cpp
# # sed -i 's|#  include <sys/wait.h>|#if defined(HAVE_SYS_WAIT_H)\n#  include <sys/wait.h>\n#endif|g' src/corelib/ncbi_process.cpp
# # sed -i 's|#  include <sys/times.h>|\/* #  include <sys/times.h> *\/|g' src/corelib/ncbi_process.cpp

# # 7. FIX: The Nuclear OS Override
# echo "Forcefully undefining Linux/Unix macros at the bottom of the config..."

# cat << 'EOF' >> "$UNIX_TEMPLATE"

# /* ========================================== */
# /* MinGW Cross-Compile Overrides              */
# /* ========================================== */
# #undef NCBI_OS_UNIX
# #undef NCBI_OS_LINUX
# #undef HOST_OS
# #undef NCBI_OS
# #define NCBI_OS "MSWIN"
# #define HOST_OS "Windows"
# EOF

# # Inject the configuration command
# echo "configure_file(\${NCBI_TREE_ROOT}/$UNIX_TEMPLATE \${NCBI_INC_ROOT}/ncbiconf_unix.h)" >> src/build-system/cmake/CMakeChecks.cmake

CUSTOM_SQLITE_LIB=$(cygpath $FOLDER_MODE $(realpath "../sqlite/libsqlite3.a"))
CUSTOM_SQLITE_INC=$(cygpath $FOLDER_MODE $(realpath "../sqlite/"))

# Helper function to check if a pattern EXISTS
verify_injected() {
    local target="$1"
    local pattern="$2"
    local fail_msg="$3"
    local grep_flags="${4:--q}" # defaults to -q (quiet)

    if ls $target 1> /dev/null 2>&1; then
        # Added -e to force grep to treat $pattern as a literal pattern, NOT a flag
        if ! grep $grep_flags -e "$pattern" $target; then
            echo -e "\n FATAL AUDIT FAILURE: $fail_msg"
            echo "   Target: $target"
            echo "   Missing Pattern: $pattern"
            exit 1
        fi
    fi
}

# Helper function to check if a pattern was DESTROYED
verify_destroyed() {
    local target="$1"
    local pattern="$2"
    local fail_msg="$3"
    local grep_flags="${4:--q}"

    if ls $target 1> /dev/null 2>&1; then
        # Added -e here as well
        if grep $grep_flags -e "$pattern" $target; then
            echo -e "\n FATAL AUDIT FAILURE: $fail_msg"
            echo "   Target: $target"
            echo "   Forbidden Pattern Found: $pattern"
            exit 1
        fi
    fi
}

UNIX_TEMPLATE=$(find src/build-system -type f \( -name "*conf*.h.in" -o -name "config*.h.in" \) | head -n 1)

if [[ ! -f "PATCHED" ]]; then
echo "Bypassing OS gates: Locating the hidden CMake config template..."

if [ -z "$UNIX_TEMPLATE" ]; then
    echo "FATAL: Could not find ANY configuration template (.h.in) in src/build-system!"
    exit 1
fi

echo "Success! Found the actual template at: $UNIX_TEMPLATE"

# echo "Forcing ncbiconf_unix.h generation for MinGW..."
# # Force CMake to generate the GCC header even on Windows
# sed -i 's/if(UNIX)/if(UNIX OR MINGW)/g' src/build-system/cmake/CMakeChecks.cmake

# Strip hardcoded MSVC flags from Windows blocks
find src/build-system/cmake -type f -name "*.cmake" -exec sed -i 's/\/MP//g; s/\/EHsc//g; s/\/Zc:__cplusplus//g' {} +

# 1. FIX: Path error (/DNDEBUG to -DNDEBUG)
echo "Patching MSVC flags for MinGW GCC compatibility..."
find src/build-system -type f -name "*.cmake" -exec sed -i 's|/DNDEBUG|-DNDEBUG|g' {} +

# 2. FIX: ssize_t clash
echo "Patching ssize_t definitions in the template..."
sed -i 's/typedef __int64 ssize_t;/\/\/ ssize_t already defined by MinGW/g' "$UNIX_TEMPLATE"

# 3. FIX: Stop CMake from defining NCBI_COMPILER_MSVC for MinGW
echo "Stripping false MSVC compiler definitions..."
sed -i 's/#define  NCBI_COMPILER_MSVC 1/\/* #undef NCBI_COMPILER_MSVC *\//g' "$UNIX_TEMPLATE"
sed -i 's/#cmakedefine NCBI_COMPILER_MSVC 1/\/* #undef NCBI_COMPILER_MSVC *\//g' "$UNIX_TEMPLATE"

# 4. FIX: Patch ncbistre.hpp so it stops trying to use Microsoft _Openprot extensions
echo "Patching ncbistre.hpp to use standard GNU libstdc++ file streams..."
sed -i 's/#if defined(NCBI_OS_MSWIN)/#if defined(NCBI_OS_MSWIN) \&\& defined(_MSC_VER)/g' include/corelib/ncbistre.hpp
sed -i 's/#ifdef NCBI_OS_MSWIN/#if defined(NCBI_OS_MSWIN) \&\& defined(_MSC_VER)/g' include/corelib/ncbistre.hpp

# 5. FIX: Destroy the hardcoded MSVC requirement gate across ALL headers
echo "Removing hardcoded MSVC restriction globally..."
find include src -type f -name "*.h" -exec sed -i 's/.*error "This toolkit is not buildable with a compiler other than MSVC.".*/\/\/ Bypassed MSVC requirement for MinGW/g' {} +

# 7. FIX: The Nuclear Header Override (OS Synchronizer)
echo "Forcefully undefining Linux macros and injecting Windows ABI at the top of the config..."
cat << 'EOF' >> "$UNIX_TEMPLATE"

/* ========================================== */
/* MinGW OS & ABI Synchronization Override    */
/* ========================================== */
#undef NCBI_OS_UNIX
#undef NCBI_OS_LINUX
#undef HOST_OS
#undef NCBI_OS
#define NCBI_OS "MSWIN"
#define HOST_OS "Windows"
#define NCBI_OS_MSWIN 1

/* Destroy MSVC flags passed by CMake command line to fix datatool */
#undef NCBI_COMPILER_MSVC

/* Disable GNU/POSIX extensions to force Windows fallbacks */
#undef HAVE_VASPRINTF
#undef HAVE_ASPRINTF
#undef HAVE_LOCALTIME_R
#undef HAVE_GMTIME_R
#undef HAVE_TIMEGM
#undef NCBI_HAVE_LOCALTIME_R
#undef NCBI_HAVE_GMTIME_R
#undef HAVE_BYTESWAP_H
#undef HAVE_ENDIAN_H
#undef HAVE_SYS_ENDIAN_H
#undef HAVE_MACHINE_ENDIAN_H
#undef HAVE_ACCEPT4
#undef HAVE_MEMRCHR
#undef HAVE_MEMCCHR
#undef HAVE_STRNDUP
#undef HAVE_SETENV
#undef HAVE_GETEUID
#undef HAVE_GETUID
#undef HAVE_SIGACTION
EOF

# 8. FIX: Fix MSVC sloppy sizeof syntax in ncbi_stack_win64.cpp
echo "Fixing non-standard sizeof syntax..."
sed -i 's/(sizeof IMAGEHLP_SYMBOL64)/(sizeof(IMAGEHLP_SYMBOL64))/g' src/corelib/ncbi_stack_win64.cpp

# 9. FIX: Disable proprietary Microsoft CRT Debug reporting
echo "Disabling MSVC-specific CRT debug hooks..."
sed -i 's/_CrtSetReportFile/\/\/ _CrtSetReportFile/g' src/corelib/ncbi_system.cpp
sed -i 's/_CrtSetReportMode/\/\/ _CrtSetReportMode/g' src/corelib/ncbi_system.cpp

# 10. FIX: Fix GCC strict function-to-data pointer casting in ncbidll.cpp
echo "Patching FARPROC to void* strict pointer casting..."
sed -i 's/entry.data = ptr;/entry.data = (void*)(intptr_t)ptr;/g' src/corelib/ncbidll.cpp

# 11. FIX: Bypass MSVC-specific fstream(FILE*) constructor in ncbifile.cpp
echo "Patching non-standard fstream(FILE*) constructor..."
sed -i 's/fstream(file)/fstream(s)/g' src/corelib/ncbifile.cpp

# 12. FIX: Remove MSVC-specific explicit std::initializer_list constructor
echo "Patching multiline non-standard std::initializer_list constructor in ncbi_pool_balancer.cpp..."
perl -0777 -pi -e 's/discrete_distribution<>\s*dd\s*\(.*?\);/discrete_distribution<> dd(weights.begin(), weights.end());/gs' src/corelib/ncbi_pool_balancer.cpp

# 13. FIX: Hijack the UNIX-specific OS file to compile the Windows implementations
echo "Hijacking ncbi_os_unix.cpp to compile Windows implementations..."
cat << 'EOF' > src/corelib/ncbi_os_unix.cpp
#include "ncbi_os_mswin.cpp"
#if __has_include("ncbi_win_security.cpp")
#include "ncbi_win_security.cpp"
#endif
EOF

# 14. FIX: Patch CityHash and FarmHash to use MinGW GCC builtins instead of Linux headers
echo "Patching 3rd-party hash libraries to use GCC builtins..."
perl -pi -e 's/#include <endian.h>/#ifndef __BYTE_ORDER\n#define __LITTLE_ENDIAN 1234\n#define __BIG_ENDIAN 4321\n#define __BYTE_ORDER __LITTLE_ENDIAN\n#endif/g' src/util/checksum/farmhash/farmhash.h

# Broadly knock out byteswap.h across the checksum directory
find src/util/checksum -type f -exec sed -i 's|#include <byteswap.h>|// #include <byteswap.h>|g' {} +

cat << 'EOF' >> src/util/checksum/cityhash/config.h
/* MinGW Cross-Compile Overrides */
#ifndef bswap_32
#define bswap_32(x) __builtin_bswap32(x)
#endif
#ifndef bswap_64
#define bswap_64(x) __builtin_bswap64(x)
#endif
EOF

# 15. FIX: Prevent mode_t collision in tar.cpp
echo "Patching mode_t collision for MinGW..."
sed -i 's/typedef unsigned int mode_t;/\/\/ typedef unsigned int mode_t;/g' src/util/compress/api/tar.cpp

# 16. FIX: Broadcast Windows OS, Unicode, and WIN32_WINNT definitions directly to all compilers via CMake
echo "Broadcasting Windows OS, Threading, and Unicode variables to compiler flags..."
cat << 'EOF' >> src/build-system/cmake/CMakeChecks.cmake

# Forcefully strip CMake's automatic POSIX thread flags
remove_definitions(-DNCBI_POSIX_THREADS -DNCBI_USE_MUTEX_PTHREADS -DNCBI_NO_THREADS)

# Globally broadcast Windows OS, Unicode, Win32 Threading, Windows 7 API, and Microsoft Threading Macro (_MT)
add_definitions(-DNCBI_OS_MSWIN=1 -DWIN32=1 -D_WIN32=1 -DUNICODE -D_UNICODE -DNCBI_WIN32_THREADS=1 -D_WIN32_WINNT=0x0601 -D_MT=1 -DNCBI_THREADS=1)
EOF

# 17. FIX: Hollow out Unix-only System V IPC and LBSM files
echo "Hollowing out Unix LBSM and IPC files..."
echo "/* Bypassed for MinGW Windows build */" > src/connect/ncbi_lbsm_ipc.c
echo "/* Bypassed for MinGW Windows build */" > src/connect/ncbi_lbsm.c
find src/connect include/connect -type f -exec sed -i 's|#include <sys/ipc.h>|// #include <sys/ipc.h>|g' {} +
find src/connect include/connect -type f -exec sed -i 's|#include <sys/shm.h>|// #include <sys/shm.h>|g' {} +

# 18. FIX: Force CMake to link Windows system libraries securely and force Console App compilation
echo "Configuring internal CMake linker for Windows libraries and Console executables..."
cat << 'EOF' >> src/build-system/cmake/CMakeChecks.cmake
# 1. Added sqlite3 to the end of this list
link_libraries(ws2_32 dbghelp advapi32 psapi iphlpapi crypt32 oleaut32 ole32 user32 netapi32 sqlite3)

set(CMAKE_CREATE_WIN32_EXE "-mconsole -mstackrealign")
set(CMAKE_CXX_CREATE_WIN32_EXE "-mconsole -mstackrealign")
set(CMAKE_C_CREATE_WIN32_EXE "-mconsole -mstackrealign")

# 2. Keep the EXE flag injection
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--allow-multiple-definition  -mconsole -mstackrealign -Wl,--subsystem,console -Wl,--demangle -Wl,--pic-executable -Wl,--support-old-code -Wl,--as-needed ")

# 3. ADD the SHARED flag injection for DLLs
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--allow-multiple-definition -Wl,--demangle -Wl,--subsystem,console -Wl,--pic-executable -Wl,-ffat-lto-objects -Wl,-flto=auto -Wl,-fstack-protector-strong -Wl,--enable-auto-import -Wl,--support-old-code -Wl,--as-needed -Wl,--dll -Wl,-shared -Wl,--no-whole-archive ")
EOF

# 19. FIX: Completely erase the MSVC-specific _Is_checked_helper block from datatool
echo "Removing MSVC STL _Is_checked_helper block from datatool..."
perl -0777 -pi -e 's/template\s*<\s*typename\s+TValue\s*,\s*typename\s+TElem\s*,\s*class\s+TTraits\s*>\s*struct\s+_Is_checked_helper.*?\{\s*\};//gs' src/serial/datatool/type.hpp

# 20. FIX: The Absolute Final Linker Boss - Overriding CMake's MinGW GUI Flags
echo "Surgically overriding CMake's internal MinGW GUI creation flags..."
cat << 'EOF' >> src/build-system/cmake/CMakeChecks.cmake

# Force CMake to translate WIN32_EXECUTABLE requests into Console apps
set(CMAKE_CXX_CREATE_WIN32_EXE "-mconsole")
set(CMAKE_C_CREATE_WIN32_EXE "-mconsole")
EOF

# 21. FIX: Adding WinMain to datatool
echo "Adding safe WinMain() proxy to datatool.cpp..."
cat << 'EOF' >> src/serial/datatool/datatool.cpp

/* MinGW Linker Hack: Safe bridge from GUI WinMain to NCBI wmain() */
#if defined(_WIN32) && defined(__MINGW32__)
#include <windows.h>
#include <shellapi.h>

extern int wmain(int argc, wchar_t** argv);

int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int nCmdShow) {
    int argc;
    // Safely fetch wide arguments directly from the Windows OS, bypassing MinGW CRT
    wchar_t** argv = CommandLineToArgvW(GetCommandLineW(), &argc);
    
    // Pass the safe array to NCBI's true entry point
    int result = wmain(argc, argv);
    
    // Clean up memory
    LocalFree(argv);
    return result;
}
#endif
EOF

# 22. FIX: GCC dllimport inline method strictness in ncbidll.hpp
echo "Moving explicitly flagged dllimport inline methods to ncbidll.cpp..."
perl -0777 -pi -e 's/const string&\s*GetName\(\)\s*const\s*\{\s*return m_Name;\s*\}/const string\& GetName() const;/g' include/corelib/ncbidll.hpp
perl -0777 -pi -e 's/const TEntries&\s*GetResolvedEntries\(\)\s*const\s*\{\s*return m_ResolvedEntries;\s*\}/const TEntries\& GetResolvedEntries() const;/g' include/corelib/ncbidll.hpp
perl -0777 -pi -e 's/TEntries&\s*GetResolvedEntries\(\)\s*\{\s*return m_ResolvedEntries;\s*\}/TEntries\& GetResolvedEntries();/g' include/corelib/ncbidll.hpp

cat << 'EOF' >> src/corelib/ncbidll.cpp

/* MinGW dllimport inline strictness fixes */
namespace ncbi {
    const string& CDll::GetName() const { return m_Name; }
    const CDllResolver::TEntries& CDllResolver::GetResolvedEntries() const { return m_ResolvedEntries; }
    CDllResolver::TEntries& CDllResolver::GetResolvedEntries() { return m_ResolvedEntries; }
}
EOF

# 23. FIX: Force CMake to use the absolute path to the freshly built datatool.exe
echo "Patching CMake to use absolute path for datatool..."
find src/build-system/cmake -type f -name "*.cmake" -exec sed -i 's/set(NCBI_DATATOOL datatool)/set(NCBI_DATATOOL "${CMAKE_BINARY_DIR}\/bin\/datatool.exe")/g' {} +

# 24. FIX: Bypass Win32 macro collision by renaming variables in parse.cpp
echo "Fixing Win32 macro collision in parse.cpp..."
sed -i 's/\bscr1\b/score1/g' src/algo/gnomon/parse.cpp
sed -i 's/\bscr2\b/score2/g' src/algo/gnomon/parse.cpp

# 25. FIX: FORCE EXPORT CGetSeqLocFromStringHelper VTABLE
echo "Hunting for CGetSeqLocFromStringHelper to patch vtable..."
# 1. Find the exact header file defining the class (searching both include/ and src/)
HELPER_HEADER=$(grep -rl "class CGetSeqLocFromStringHelper" src/ include/ | grep "\.hpp$" | head -n 1)
if [ ! -z "$HELPER_HEADER" ]; then
    echo "Found target header: $HELPER_HEADER"
    # 2. Use sed to inject __declspec(dllexport) immediately after the 'class' keyword
    # This turns: class CGetSeqLocFromStringHelper
    # Into:       class __declspec(dllexport) CGetSeqLocFromStringHelper
    sed -i 's/class CGetSeqLocFromStringHelper/class __declspec(dllexport) CGetSeqLocFromStringHelper/g' "$HELPER_HEADER"
    echo "Successfully patched $HELPER_HEADER!"
else
    echo "WARNING: Could not find CGetSeqLocFromStringHelper declaration."
fi

# 26. FIX: Forcefully override CMake's find_package(SQLite3) behavior
echo "Forcefully overriding CMake's find_package(SQLite3) behavior..."
cat << EOF >> src/build-system/cmake/CMakeChecks.cmake
set(SQLITE3_LIBRARY "${CUSTOM_SQLITE_LIB}" CACHE FILEPATH "Forced SQLite" FORCE)
set(SQLite3_LIBRARY "${CUSTOM_SQLITE_LIB}" CACHE FILEPATH "Forced SQLite" FORCE)
set(SQLITE3_INCLUDE_DIR "${CUSTOM_SQLITE_INC}" CACHE PATH "Forced SQLite Inc" FORCE)
set(SQLite3_INCLUDE_DIR "${CUSTOM_SQLITE_INC}" CACHE PATH "Forced SQLite Inc" FORCE)
EOF

# 27. FIX: Patch mutex and vtable DLL exports for some headers and sources
echo "Hunting for GenBank loader locks to patch mutex and vtable exports..."
# Search inside both src/ and include/ for the base lock class declaration
GBL_LOCK_FILE=$(grep -rl "class CInfoLock_Base" src/ include/ | head -n 1)
if [ ! -z "$GBL_LOCK_FILE" ]; then
    echo "Found GenBank lock file: $GBL_LOCK_FILE"
    # 1. Patch CInfoLock_Base
    sed -i 's/class NCBI_[A-Z_]*_EXPORT CInfoLock_Base/class CInfoLock_Base/g' "$GBL_LOCK_FILE"
    sed -i 's/class CInfoLock_Base/class __declspec(dllexport) CInfoLock_Base/g' "$GBL_LOCK_FILE"
    # 2. Patch CInfoRequestorLock (usually in the same file)
    sed -i 's/class NCBI_[A-Z_]*_EXPORT CInfoRequestorLock/class CInfoRequestorLock/g' "$GBL_LOCK_FILE"
    sed -i 's/class CInfoRequestorLock/class __declspec(dllexport) CInfoRequestorLock/g' "$GBL_LOCK_FILE"  
    echo "Successfully patched GenBank locks in $GBL_LOCK_FILE!"
else
    echo "WARNING: Could not find GenBank lock declarations."
fi

# 28. FIX: Patch FreeTDS broken headers 
echo "Patching FreeTDS broken headers for MinGW..."
# 1. Fix the getpid() macro conflict
SYSDEP_FILE="include/dbapi/driver/ftds14/freetds/freetds/sysdep_private.h"
if [ -f "$SYSDEP_FILE" ]; then
    # Delete the line containing the broken getpid macro
    sed -i '/#define getpid()/d' "$SYSDEP_FILE"
    echo "Successfully patched $SYSDEP_FILE!"
fi

# 29. Fix the missing Unix readpassphrase header
READPASS_FILE="include/dbapi/driver/ftds14/freetds/freetds/replacements/readpassphrase.h"
if [ -f "$READPASS_FILE" ]; then
    # Delete the line attempting to include the non-existent system header
    sed -i '/# include <readpassphrase.h>/d' "$READPASS_FILE"
    echo "Successfully patched $READPASS_FILE!"
fi

# 30. Fix the DBAPI macro collision
echo "Patching DBAPI Windows interface macro collision..."
DBAPI_UTILS_FILE="include/dbapi/driver/impl/dbapi_driver_utils.hpp"
if [ -f "$DBAPI_UTILS_FILE" ]; then
    # Inject #undef interface at the very beginning of the file using sed
    sed -i '1s/^/#ifdef interface\n#undef interface\n#endif\n\n/' "$DBAPI_UTILS_FILE"
    echo "Successfully patched $DBAPI_UTILS_FILE!"
fi

# 31. Fix the missing strlcpy for FreeTDS on Windows
echo "Fixing missing strlcpy for FreeTDS on Windows..."
cat << 'EOF' >> src/build-system/cmake/CMakeChecks.cmake
# Dynamically replace strlcpy with snprintf to guarantee null-termination on Windows
add_compile_options("-Dstrlcpy(d,s,l)=snprintf(d,l,\"%s\",s)")
EOF

# 32. Fix the GCC 14 strict pointer errors for FreeTDS
echo "Fixing GCC 14 strict pointer errors for FreeTDS..."
cat << 'EOF' >> src/build-system/cmake/CMakeChecks.cmake
# Downgrade pointer mismatches from Fatal Errors back to Warnings
add_compile_options("-Wno-error=incompatible-pointer-types")
EOF

# 33. Fix the 'interface' macro across ALL DBAPI headers
echo "Nuking 'interface' macro across ALL DBAPI headers..."
# Find all headers in the impl folder and rename the variable globally
find src/dbapi/driver -type f \( -name "*.cpp" -o -name "*.hpp" \) -exec sed -i 's/\binterface\b/iface/g' {} +
find include/dbapi/driver -type f \( -name "*.cpp" -o -name "*.hpp" \) -exec sed -i 's/\binterface\b/iface/g' {} +

# 34. Fix the FreeTDS Unix networking headers on Windows
echo "Fixing FreeTDS Unix networking headers on Windows..."
# Comment out the Unix-specific netinet/tcp.h include
sed -i 's|#include <netinet/tcp.h>|/* #include <netinet/tcp.h> */|g' src/dbapi/driver/ftds14/freetds/tds/net.c

# Prepend Windows-specific macro fallbacks to the top of net.c
cat << 'EOF' > patch_net.c
#ifdef _WIN32
#include <winsock2.h>
#ifndef poll
#define poll WSAPoll
#endif
#define socketpair(d,t,p,s) (-1)
#ifndef AF_UNIX
#define AF_UNIX 1
#endif
#endif
EOF

# Stitch our patch and the original file together, then replace the original
cat src/dbapi/driver/ftds14/freetds/tds/net.c >> patch_net.c
mv patch_net.c src/dbapi/driver/ftds14/freetds/tds/net.c

# 35. Fix the FreeTDS obsolete condition variable reassignments
echo "Fixing FreeTDS obsolete condition variable reassignments..."
# Comment out the illegal function pointer assignments by prepending //
sed -i 's/tds_raw_cond_init[ \t]*=/\/\/ &/' src/dbapi/driver/ftds14/freetds/utils/tds_cond.c
sed -i 's/tds_raw_cond_destroy[ \t]*=/\/\/ &/' src/dbapi/driver/ftds14/freetds/utils/tds_cond.c
sed -i 's/tds_raw_cond_signal[ \t]*=/\/\/ &/' src/dbapi/driver/ftds14/freetds/utils/tds_cond.c
sed -i 's/tds_raw_cond_timedwait[ \t]*=/\/\/ &/' src/dbapi/driver/ftds14/freetds/utils/tds_cond.c

# 36. Fix the missing Windows headers in FreeTDS utilities
echo "Injecting missing Windows headers into FreeTDS utilities..."
sed -i '1i #include <windows.h>' src/dbapi/driver/ftds14/freetds/utils/tds_cond.c
sed -i '1i #include <windows.h>' src/dbapi/driver/ftds14/freetds/utils/win_mutex.c
sed -i '1i #include <windows.h>' src/dbapi/driver/ftds14/freetds/utils/smp.c

# 37. Fix the FreeTDS function pointer vs. macro collisions in tds_cond.c
echo "Fixing FreeTDS function pointer vs. macro collisions in tds_cond.c..."
sed -i '/int (\*tds_raw_cond_init)/c\int tds_raw_cond_init(tds_condition * cond) { return new_cond_init(cond); }' src/dbapi/driver/ftds14/freetds/utils/tds_cond.c
sed -i '/int (\*tds_raw_cond_destroy)/c\int tds_raw_cond_destroy(tds_condition * cond) { return new_cond_destroy(cond); }' src/dbapi/driver/ftds14/freetds/utils/tds_cond.c
sed -i '/int (\*tds_raw_cond_signal)/c\int tds_raw_cond_signal(tds_condition * cond) { return new_cond_signal(cond); }' src/dbapi/driver/ftds14/freetds/utils/tds_cond.c
sed -i '/int (\*tds_raw_cond_timedwait)/c\int tds_raw_cond_timedwait(tds_condition * cond, tds_raw_mutex * mtx, int timeout_sec) { return new_cond_timedwait(cond, mtx, timeout_sec); }' src/dbapi/driver/ftds14/freetds/utils/tds_cond.c

# 38. Fix the win_mutex.c inline redefinition conflict
echo "Resolving win_mutex.c inline redefinition conflict..."
sed -i '/tds_raw_mutex_trylock(tds_raw_mutex/c\int tds_raw_mutex_trylock_IGNORED(tds_raw_mutex * mutex)' src/dbapi/driver/ftds14/freetds/utils/win_mutex.c

# 39. FIX: Resolve strict pointer casting errors in Windows API hooking
echo "Patching explicit pointer casts for MinGW in ncbi_win_hook.cpp..."
# 1. Cast internal hook addresses to const void* and LoadLibraryA to PVOID to satisfy GCC's stricter type checking
sed -i -e 's/IsPatched(sm_/IsPatched((const void*)sm_/g' \
       -e 's/ModuleFromAddress(CKernell32::LoadLibraryA)/ModuleFromAddress((PVOID)CKernell32::LoadLibraryA)/g' \
       src/dbapi/driver/ncbi_win_hook.cpp

# 40. FIX: Align FreeTDS condition variables with Native Windows API
echo "Patching FreeTDS condition variables for Windows compatibility..."
# 1. Map FreeTDS condition variables directly to native Windows CONDITION_VARIABLE and fix macro collisions
sed -i -e 's/TDS_CONDITION_VARIABLE/CONDITION_VARIABLE/g' \
       -e 's/int tds_raw_cond_destroy(tds_condition/#undef tds_raw_cond_destroy\nint tds_raw_cond_destroy(tds_condition/g' \
       -e 's/int tds_raw_cond_signal(tds_condition/#undef tds_raw_cond_signal\nint tds_raw_cond_signal(tds_condition/g' \
       src/dbapi/driver/ftds14/freetds/utils/tds_cond.c

# 41. FIX: Enforce minimum Windows API version for FreeTDS
echo "Forcing _WIN32_WINNT to 0x0600 (Windows Vista) in FreeTDS windows.h..."
# 1. Undefine any existing OS targets and lock the API level to ensure required networking/threading features are exposed
sed -i '1s/^/#undef _WIN32_WINNT\n#define _WIN32_WINNT 0x0600\n/' include/dbapi/driver/ftds14/freetds/freetds/windows.h              

# 42. FIX: Resolve variable name collisions in GenBank flatfile parser
echo "Renaming conflicting variables in genref.cpp..."
# 1. Change 'grp1' and 'grp2' to 'gref1' and 'gref2' to avoid shadowing/clashes with MinGW system macros
sed -i -e 's/\bgrp1\b/gref1/g' -e 's/\bgrp2\b/gref2/g' src/objtools/flatfile/genref.cpp

# 43. FIX: Ensure Windows Sockets API is loaded early for FreeTDS
echo "Injecting winsock2.h into FreeTDS windows.h..."
# 1. Prepend the winsock2 header to provide necessary networking definitions before FreeTDS types are declared
sed -i '1i #include <winsock2.h>' include/dbapi/driver/ftds14/freetds/freetds/windows.h

# 44. FIX: Disable POSIX threading in FreeTDS on Windows
echo "Stripping POSIX pthread macros from FreeTDS headers..."
# 1. Find all FreeTDS headers and comment out TDS_HAVE_PTHREAD to force native Windows threading fallbacks
find include/dbapi/driver/ftds14 src/dbapi/driver/ftds14 -type f -name "*.h" -exec sed -i 's/#define TDS_HAVE_PTHREAD/\/\/ #define TDS_HAVE_PTHREAD/g' {} +

# 45. FIX: Nuke redundant FreeTDS Windows mutex implementations
echo "Disabling compilation of FreeTDS win_mutex.c..."
# 1. Wrap the entire file in an #if 0 block to avoid multiple definition linker errors against standard libraries
sed -i '1s/^/#if 0\n/' src/dbapi/driver/ftds14/freetds/utils/win_mutex.c
echo "#endif" >> src/dbapi/driver/ftds14/freetds/utils/win_mutex.c

# 46. FIX: Nuke redundant FreeTDS condition variable implementations
echo "Disabling compilation of FreeTDS tds_cond.c..."
# 1. Wrap the entire file in an #if 0 block (superseding previous patches to guarantee it doesn't conflict during linking)
sed -i '1s/^/#if 0\n/' src/dbapi/driver/ftds14/freetds/utils/tds_cond.c
echo "#endif" >> src/dbapi/driver/ftds14/freetds/utils/tds_cond.c

# 47. FIX: MinGW GCC Nested Class Export Bug for CCachedRowInfo::SInfo
echo "Patching CCachedRowInfo::SInfo for explicit DLL export..."
# 1. Target the exact file
TARGET_FILE="include/dbapi/driver/impl/dbapi_driver_utils.hpp"
# 2. Dynamically grab the exact export macro used in this header (usually NCBI_DBAPIDRIVER_EXPORT)
EXPORT_MACRO=$(grep -o -m 1 "NCBI_[A-Z_]*_EXPORT" "$TARGET_FILE")
echo "Found export macro: $EXPORT_MACRO"
# 3. Patch the nested struct to ensure GCC exports the constructor to the DLL
sed -i "s/struct SInfo/struct $EXPORT_MACRO SInfo/g" "$TARGET_FILE"
# # 4. Verify the patch was applied
# grep -n "struct $EXPORT_MACRO SInfo" "$TARGET_FILE"

# 48. FIX: Update makeblastdb dependencies to build static libs
echo "Updating makeblastdb dependencies to fix static linking"
sed -i 's/NCBI_uses_toolkit_libraries(blastinput writedb)/NCBI_uses_toolkit_libraries(blastinput writedb blastdb blast xobjmgr general xser xutil xncbi)/g' src/app/blastdb/CMakeLists.makeblastdb.app.txt

# 49. FIX: Inject the configuration command
echo "Injecting configuration command to force ncbiconf_unix.h generation"
echo "configure_file(\${NCBI_TREE_ROOT}/$UNIX_TEMPLATE \${NCBI_INC_ROOT}/ncbiconf_unix.h)" >> src/build-system/cmake/CMakeChecks.cmake

touch PATCHED
fi

echo "==================================================="
echo "AUDIT PHASE: Verifying OS Gates and ABI Injections..."
echo "==================================================="
# --- 1. CMake System Checks ---
# verify_injected "src/build-system/cmake/CMakeChecks.cmake" "if(UNIX OR MINGW)" "CMake UNIX gate bypass failed."
verify_destroyed "src/build-system/cmake/*.cmake" "/DNDEBUG" "MSVC /DNDEBUG flags were not successfully stripped." -rq
verify_injected "src/build-system/cmake/CMakeChecks.cmake" "remove_definitions(-DNCBI_POSIX_THREADS" "POSIX thread strip failed."
verify_injected "src/build-system/cmake/CMakeChecks.cmake" "-DNCBI_WIN32_THREADS=1" "Windows OS variables broadcast failed."
verify_injected "src/build-system/cmake/CMakeChecks.cmake" "sqlite3" "CMake linker modifications missing SQLite3."
verify_injected "src/build-system/cmake/CMakeChecks.cmake" "set(CMAKE_CXX_CREATE_WIN32_EXE \"-mconsole\")" "GUI override to Console failed."
verify_injected "src/build-system/cmake/CMakeChecks.cmake" "-Wno-error=incompatible-pointer-types" "GCC 14 strict pointer override failed."
# --- 2. Unix Template & OS Header Checks ---
verify_injected "$UNIX_TEMPLATE" "ssize_t already defined" "ssize_t collision not patched in UNIX template."
verify_injected "$UNIX_TEMPLATE" "MinGW OS & ABI Synchronization Override" "Nuclear header override failed."
verify_injected "include/corelib/ncbistre.hpp" "defined(NCBI_OS_MSWIN) && defined(_MSC_VER)" "ncbistre.hpp Microsoft extensions not disabled."
# verify_destroyed "include/ src/" "This toolkit is not buildable with a compiler other than MSVC" "Hardcoded MSVC requirement still exists!" -rq
verify_injected "src/corelib/ncbi_os_unix.cpp" "ncbi_os_mswin.cpp" "UNIX OS file hijack failed."
# --- 3. Corelib & Checksum Logic Checks ---
verify_injected "src/corelib/ncbi_stack_win64.cpp" "sizeof(IMAGEHLP_SYMBOL64)" "sizeof syntax not fixed."
verify_injected "src/corelib/ncbidll.cpp" "(void\*)(intptr_t)ptr;" "FARPROC pointer casting not patched."
verify_injected "src/corelib/ncbifile.cpp" "fstream(s)" "fstream(FILE*) bypass failed."
verify_injected "src/corelib/ncbi_pool_balancer.cpp" "dd(weights.begin(), weights.end())" "std::initializer_list constructor bypass failed."
verify_injected "src/util/checksum/cityhash/config.h" "__builtin_bswap" "CityHash GCC builtin overrides failed."
verify_injected "src/util/compress/api/tar.cpp" "// typedef unsigned int mode_t;" "mode_t collision bypass failed."
verify_injected "include/corelib/ncbidll.hpp" "const string& GetName() const;" "dllimport inline strictness bypass failed."
# --- 4. Datatool & Serial Checks ---
verify_injected "src/serial/datatool/datatool.cpp" "WINAPI WinMain" "Safe wmain bridge missing in datatool."
verify_destroyed "src/serial/datatool/type.hpp" "struct _Is_checked_helper" "MSVC STL checked helper block not deleted."
verify_injected "src/algo/gnomon/parse.cpp" "\bscore1\b" "Win32 macro collision (scr1/scr2) not patched." -Eq
# --- 5. DLL Export VTable Checks ---
verify_injected "include/ src/" "class __declspec(dllexport) CGetSeqLocFromStringHelper" "CGetSeqLocFromStringHelper vtable export patch failed." -rq
verify_injected "include/ src/" "class __declspec(dllexport) CInfoLock_Base" "GenBank locks vtable export patch failed." -rq
# --- 6. FreeTDS Dependency Checks ---
verify_destroyed "include/dbapi/driver/ftds14/freetds/freetds/sysdep_private.h" "#define getpid()" "FreeTDS getpid() macro collision not resolved."
# verify_injected "include/dbapi/driver/impl/dbapi_driver_utils.hpp" "#undef interface" "DBAPI 'interface' collision undefine failed."
verify_destroyed "src/dbapi/driver/ include/dbapi/driver/" "^[[:space:]]*#define TDS_HAVE_PTHREAD" "POSIX pthread macros not stripped from FreeTDS." -rqE
verify_injected "src/dbapi/driver/ftds14/freetds/tds/net.c" "#ifdef _WIN32" "FreeTDS Unix networking fallback not prepended."
verify_injected "src/dbapi/driver/ftds14/freetds/utils/tds_cond.c" "CONDITION_VARIABLE" "FreeTDS condition variables not aligned with Windows API."
verify_injected "include/dbapi/driver/ftds14/freetds/freetds/windows.h" "#define _WIN32_WINNT 0x0600" "FreeTDS Win32 API level override failed."
verify_injected "src/dbapi/driver/ftds14/freetds/utils/win_mutex.c" "^#if 0" "FreeTDS win_mutex.c not disabled." -Eq
# --- 7. Final Configuration Output Checks ---
verify_injected "src/app/blastdb/CMakeLists.makeblastdb.app.txt" "xncbi" "makeblastdb static linker dependencies not injected."
verify_injected "src/build-system/cmake/CMakeChecks.cmake" "configure_file.*ncbiconf_unix.h" "UNIX config generation override missing." -Eq
echo " All OS Gates bypassed successfully. Patches strictly verified."
echo "==================================================="

rm -f src/build-system/cmake/toolchains/$toolchain_file

touch src/build-system/cmake/toolchains/$toolchain_file

#############GCC
echo "set(NCBI_PTBCFG_FLAGS_DEFINED YES)" >>src/build-system/cmake/toolchains/$toolchain_file
echo "include_guard(GLOBAL)" >>src/build-system/cmake/toolchains/$toolchain_file
#-gdwarf-4
#-Wall
#-std=gnu2x
#-std=gnu++20
#-O3 -fPIC
echo "set(CMAKE_C_FLAGS_INIT         \"-Wno-format-y2k -Wno-date-time -Wno-attributes -Wno-unused-parameter -Wno-ignored-attributes -Wa,-mbig-obj -pedantic $SHLIB_OPENMP_CXXFLAGS -pthread -ffat-lto-objects -flto=auto -fstack-protector-strong -static-libgcc -march=native -mtune=generic -maccumulate-outgoing-args -mconsole -mstackrealign \")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(CMAKE_CXX_FLAGS_INIT       \"-Wno-format-y2k -Wno-date-time -Wno-attributes -Wno-unused-parameter -Wno-ignored-attributes -Wa,-mbig-obj -DLIBICONV_STATIC -pedantic $SHLIB_OPENMP_CXXFLAGS -pthread -ffat-lto-objects -flto=auto -fstack-protector-strong -static-libstdc++ -march=native -mtune=generic -maccumulate-outgoing-args -mconsole -mstackrealign \")" >>src/build-system/cmake/toolchains/$toolchain_file
# echo "set(NCBI_COMPILER_FLAGS_SSE       \"-msse4.2\")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(NCBI_COMPILER_FLAGS_SSE       \"-march=native\")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(NCBI_COMPILER_FLAGS_COVERAGE  \"--coverage\")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(NCBI_LINKER_FLAGS_COVERAGE     \"--coverage\")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(NCBI_COMPILER_FLAGS_MAXDEBUG  \"-fsanitize=address -fstack-check\")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(NCBI_LINKER_FLAGS_MAXDEBUG   \"-fsanitize=address\")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(NCBI_LINKER_FLAGS_STATICCOMPONENTS \"-static-libgcc -static-libstdc++\")" >>src/build-system/cmake/toolchains/$toolchain_file

echo "set(CMAKE_HOST_SYSTEM_NAME Windows)" >>src/build-system/cmake/toolchains/$toolchain_file
# # echo "set(CMAKE_HOST_SYSTEM_NAME mingw)" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(CMAKE_SYSTEM_NAME Linux)" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(CMAKE_BUILD_TYPE Release)" >>src/build-system/cmake/toolchains/$toolchain_file
# echo "set(BUILD_SHARED_LIBS ON)" >>src/build-system/cmake/toolchains/$toolchain_file
 ##echo "set(NCBI_DLL_BUILD ON)" >>src/build-system/cmake/toolchains/$toolchain_file
 echo "set(BUILD_SHARED_LIBS OFF)" >>src/build-system/cmake/toolchains/$toolchain_file
 echo "set(BUILD_STATIC_LIBS ON)" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(CMAKE_POSITION_INDEPENDENT_CODE ON)" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(NCBI_COMPILER GNU)" >>src/build-system/cmake/toolchains/$toolchain_file

echo "set(NCBI_COMPILER_VERSION $COMPILER_VER_NOPUNCT)" >>src/build-system/cmake/toolchains/$toolchain_file

echo "set(CMAKE_STATIC_LIBRARY_SUFFIX_C \".a\")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(CMAKE_STATIC_LIBRARY_PREFIX_C \"lib\")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(CMAKE_STATIC_LIBRARY_SUFFIX_CXX \".a\")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(CMAKE_STATIC_LIBRARY_PREFIX_CXX \"lib\")" >>src/build-system/cmake/toolchains/$toolchain_file

echo "set(CMAKE_SHARED_LIBRARY_SUFFIX_C \".dll\")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(CMAKE_SHARED_LIBRARY_PREFIX_C \"lib\")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(CMAKE_SHARED_LIBRARY_SUFFIX_CXX \".dll\")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(CMAKE_SHARED_LIBRARY_PREFIX_CXX \"lib\")" >>src/build-system/cmake/toolchains/$toolchain_file

echo "set(CMAKE_SYSROOT_LINK \"/x86_64-w64-mingw32.static.posix/lib\")" >> src/build-system/cmake/toolchains/$toolchain_file
echo "set(CMAKE_SYSROOT_COMPILE \"/x86_64-w64-mingw32.static.posix/include\")" >> src/build-system/cmake/toolchains/$toolchain_file
echo "set(CMAKE_SYSROOT \"/x86_64-w64-mingw32.static.posix\")" >> src/build-system/cmake/toolchains/$toolchain_file

# # echo "set(CMAKE_C_COMPILER \"$CC\" CACHE FILEPATH \"Force C\" FORCE)" >> src/build-system/cmake/toolchains/$toolchain_file
# # echo "set(CMAKE_CXX_COMPILER \"$CXX\" CACHE FILEPATH \"Force C\" FORCE)" >> src/build-system/cmake/toolchains/$toolchain_file
# # #echo "set(CMAKE_C_COMPILER_LINKER \"$C_LD\" CACHE FILEPATH \"Force C\" FORCE)" >> src/build-system/cmake/toolchains/$toolchain_file
# # #echo "set(CMAKE_CXX_COMPILER_LINKER \"$CXX_LD\" CACHE FILEPATH \"Force C\" FORCE)" >> src/build-system/cmake/toolchains/$toolchain_file
# # #echo "set(LD \"$CXX_LD\")" >> src/build-system/cmake/toolchains/$toolchain_file
# # 2. Force CMake to use these exact forward-slashed tools so it NEVER scrapes GCC's backslashes
# echo "set(CMAKE_C_COMPILER \"$CC\" CACHE FILEPATH \"Force C\" FORCE)" >> src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_CXX_COMPILER \"$CXX\" CACHE FILEPATH \"Force CXX\" FORCE)" >> src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_LINKER \"$C_LD\" CACHE FILEPATH \"Force Linker\" FORCE)" >> src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_C_COMPILER_LINKER \"$C_LD\" CACHE FILEPATH \"Force C Linker\" FORCE)" >> src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_CXX_COMPILER_LINKER \"$CXX_LD\" CACHE FILEPATH \"Force CXX Linker\" FORCE)" >> src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_AR \"$AR_PATH\" CACHE FILEPATH \"Force AR\" FORCE)" >> src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_NM \"$NM_PATH\" CACHE FILEPATH \"Force NM\" FORCE)" >> src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_RANLIB \"$RANLIB_PATH\" CACHE FILEPATH \"Force RANLIB\" FORCE)" >> src/build-system/cmake/toolchains/$toolchain_file



#-U_WIN32
#-ftrack-macro-expansion=0 -pipe 
#-std=gnu2x
#-O3 -fPIC
#############COMPILER - MSYS2 - GCC#################################
echo "set(CMAKE_C_FLAGS_RELEASE \" -I'$NCBI_DIR/BUILD/inc/common/config' -Wformat -Werror=format-security -DHAVE_IOSTREAM=1 -DNCBI_OS_OSF1=1 -D_WIN32 -D_WIN64 -DHAVE_INTTYPES_H=1 -DHAVE_NETINET_TCP_H=1 -Wa,-mbig-obj -D_FORTIFY_SOURCE=2 -pedantic $SHLIB_OPENMP_CXXFLAGS -pthread -static-libgcc -march=native -mtune=generic -maccumulate-outgoing-args -Wno-format-y2k -Wno-date-time -Wno-attributes -mconsole -mstackrealign \")" >>src/build-system/cmake/toolchains/$toolchain_file
#-Wall
#-D_DEBUG
#-UNDEBUG
#-fstack-protector-strong
#-fdiagnostics-color=always
#-U_WIN32
#-Wall
#-ftrack-macro-expansion=0 -pipe 
#-std=gnu++20
#-O3 -fPIC
echo "set(CMAKE_CXX_FLAGS_RELEASE \" -Wformat -Werror=format-security -Wdate-time -DHAVE_IOSTREAM=1 -DNCBI_OS_OSF1=1 -D_WIN32 -D_WIN64 -DHAVE_INTTYPES_H=1 -DHAVE_NETINET_TCP_H=1 -D_FORTIFY_SOURCE=2 -DLIBICONV_STATIC -pedantic $SHLIB_OPENMP_CXXFLAGS -pthread -ffat-lto-objects -flto=auto -fstack-protector-strong -static-libstdc++ -march=native -mtune=generic -maccumulate-outgoing-args -Wno-format-y2k -Wno-date-time -Wno-attributes -Wa,-mbig-obj -mconsole -mstackrealign \")" >>src/build-system/cmake/toolchains/$toolchain_file
###################################################

##############LINKER - GCC##############
echo "set(CMAKE_EXE_LINKER_FLAGS_INIT  \" \${CMAKE_EXE_LINKER_FLAGS} -Wl,--allow-multiple-definition -Wl,--demangle -mconsole -mstackrealign -Wl,--subsystem,console -Wl,--pic-executable -Wl,--support-old-code -Wl,--as-needed -Wl,--start-group -lws2_32 -ldbghelp -ladvapi32 -lpsapi -liphlpapi -lcrypt32  ${CUSTOM_SQLITE_LIB} -liconv -lintl -lbcrypt -Wl,--end-group \")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(CMAKE_SHARED_LINKER_FLAGS_INIT  \" \${CMAKE_SHARED_LINKER_FLAGS} -Wl,--allow-multiple-definition -Wl,--demangle -Wl,--subsystem,console -Wl,--pic-executable -Wl,-ffat-lto-objects -Wl,-flto=auto -Wl,-fstack-protector-strong -Wl,--enable-auto-import -Wl,--support-old-code -Wl,--as-needed -Wl,--dll -Wl,-shared -Wl,--no-whole-archive -Wl,--start-group -lws2_32 -ldbghelp -ladvapi32 -lpsapi -liphlpapi -lcrypt32 ${CUSTOM_SQLITE_LIB} -liconv -lintl -lbcrypt -Wl,--end-group \")" >>src/build-system/cmake/toolchains/$toolchain_file

echo "set(CMAKE_BUILD_WITH_INSTALL_RPATH ON)" >>src/build-system/cmake/toolchains/$toolchain_file
#echo "set(CMAKE_INSTALL_RPATH='\\$$ORIGIN'\")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(CMAKE_INSTALL_RPATH_USE_LINK_PATH OFF)" >>src/build-system/cmake/toolchains/$toolchain_file

# echo "set(CMAKE_LINKER_TYPE LLD)" >>src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_C_USING_LINKER_LLD \"ld.exe\")" >>src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_C_USING_LINKER_MODE TOOL)" >>src/build-system/cmake/toolchains/$toolchain_file
################################

echo "include(\"$NCBI_DIR/ncbi_components.cmake\")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(CMAKE_CXX_EXTENSIONS ON)" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(MINGW32 1)" >>src/build-system/cmake/toolchains/$toolchain_file

echo "set(NCBI_COMPILER_VERSION_DOTTED \"$COMPILER_VER\")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(NCBI_COMPILER_VERSION \"$COMPILER_VER_NOPUNCT\")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(MSVC_VERSION \"$COMPILER_VER_NOPUNCT\")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(_MSC_VER \"$COMPILER_VER_NOPUNCT\")" >>src/build-system/cmake/toolchains/$toolchain_file

# echo "Patching NCBI cURL array bug..."
# sed -i 's/\${NCBI_COMPONENT_CURL_LIBS}/"${CURL_LIBRARY}"/g' src/build-system/cmake/CMake.NCBIComponentsUNIXex.cmake
# sed -i 's/\${CURL_LIBRARIES}/"${CURL_LIBRARY}"/g' src/build-system/cmake/CMake.NCBIComponentsUNIXex.cmake

# echo 'set(WIN32 0)' >> src/build-system/cmake/toolchains/$toolchain_file
# echo 'set(UNIX 1)' >> src/build-system/cmake/toolchains/$toolchain_file
# echo 'set(CYGWIN 0)' >> src/build-system/cmake/toolchains/$toolchain_file
# echo 'set(MSYS 1)' >> src/build-system/cmake/toolchains/$toolchain_file
# echo 'set(HOST_OS windows)' >> src/build-system/cmake/toolchains/$toolchain_file

echo 'set(WIN32 0)' >> src/build-system/cmake/toolchains/$toolchain_file
echo 'set(UNIX 1)' >> src/build-system/cmake/toolchains/$toolchain_file
echo 'set(CYGWIN 0)' >> src/build-system/cmake/toolchains/$toolchain_file
echo 'set(MSYS 0)' >> src/build-system/cmake/toolchains/$toolchain_file
echo 'set(CMAKE_SYSTEM_VERSION "1.0.0")' >> src/build-system/cmake/toolchains/$toolchain_file
echo 'set(CMAKE_HOST_SYSTEM_VERSION "1.0.0")' >> src/build-system/cmake/toolchains/$toolchain_file

# --- ADD THIS LINE TO FIX TRY_COMPILE ---
echo 'set(CMAKE_EXECUTABLE_SUFFIX ".exe")' >> src/build-system/cmake/toolchains/$toolchain_file
# ----------------------------------------

echo 'set($ENV{OSTYPE} linux)' >> src/build-system/cmake/toolchains/$toolchain_file
echo 'set(HOST_OS Linux)' >> src/build-system/cmake/toolchains/$toolchain_file

# echo "message('CC: \${CC}')" >>src/build-system/cmake/toolchains/$toolchain_file
# echo "message('CXX: \${CXX}')" >>src/build-system/cmake/toolchains/$toolchain_file
# echo "message('CMAKE_C_COMPILER: \${CMAKE_C_COMPILER}')" >>src/build-system/cmake/toolchains/$toolchain_file
# echo "message('CMAKE_CXX_COMPILER: \${CMAKE_CXX_COMPILER}')" >>src/build-system/cmake/toolchains/$toolchain_file
# echo "message('CMAKE_C_COMPILER_LINKER: \${CMAKE_C_COMPILER_LINKER}')" >>src/build-system/cmake/toolchains/$toolchain_file
# echo "message('CMAKE_CXX_COMPILER_LINKER: \${CMAKE_CXX_COMPILER_LINKER}')" >>src/build-system/cmake/toolchains/$toolchain_file

##################################################################################################

if [[ ! -s $SCRIPTPATH/../sqlite/libsqlite3.a ]]; then
# echo "Fetching latest SQLite version info..."

# # 1. Download the HTML page silently (-q) and output it to a file (-O)
# wget -qO download.html https://sqlite.org/download.html

# # 2. Extract the relative path (e.g., "2026/sqlite-amalgamation-3510300.zip")
# SQLITE_REL_PATH=$(grep -Eo "[0-9]{4}/sqlite-amalgamation-[0-9]+\.zip" download.html | head -n 1)

# if [ -z "$SQLITE_REL_PATH" ]; then
#     echo "Error: Could not parse SQLite download path."
#     exit 1
# fi

# echo "Latest SQLite path found: $SQLITE_REL_PATH"

# # 3. Download the actual amalgamation zip
# wget "https://sqlite.org/${SQLITE_REL_PATH}"
cd ../sqlite/
# # 4. Extract and clean up
unzip ../sqlite/sqlite.zip
rm -f download.html sqlite.zip

cd sqlite-amalgamation-*

export PATH="$RTOOLS_BIN:$PATH"

#-std=gnu2x -O3 -fPIC

gcc $PKG_CFLAGS $SHLIB_OPENMP_CXXFLAGS -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -pedantic $SHLIB_OPENMP_CXXFLAGS -pthread -static-libgcc -march=native -mtune=generic -maccumulate-outgoing-args -Wno-format-y2k -mconsole -mstackrealign  \
  -DSQLITE_ENABLE_UNLOCK_NOTIFY=1 \
  -DSQLITE_ENABLE_MEMSYS5 \
  -DSQLITE_ENABLE_MEMSYS3 \
  -DSQLITE_ENABLE_FTS3 \
  -DSQLITE_ENABLE_FTS4 \
  -DSQLITE_ENABLE_FTS5 \
  -DSQLITE_ENABLE_UPDATE_DELETE_LIMIT \
  -DSQLITE_ENABLE_COLUMN_METADATA \
  -DSQLITE_ENABLE_RTREE \
  -DSQLITE_ENABLE_GEOPOLY \
  -DSQLITE_ENABLE_SESSION \
  -DSQLITE_ENABLE_PREUPDATE_HOOK \
  -DSQLITE_ENABLE_DBPAGE_VTAB \
  -DSQLITE_ENABLE_DBSTAT_VTAB \
  -c sqlite3.c -o sqlite3.o

# 3. Archive it into a static library
ar rcs libsqlite3.a sqlite3.o

# mkdir -p ../sqlite/
cp -f libsqlite3.a ../
fi