#!/bin/bash

#pacman -Syuu --noconfirm

#Get script directory from mingw shell for sysroot
SCRIPTPATH="$(
    cd -- "$(dirname "$0")" >/dev/null 2>&1
    pwd -P
)"
## NCBI_DIR=`pwd`
# NCBI_DIR=$(cygpath -u $(pwd))
FOLDER_MODE="-u" #for cygpath
NCBI_DIR=$(cygpath $FOLDER_MODE $(pwd))
SCRIPTPATH=$(cygpath $FOLDER_MODE $SCRIPTPATH)
DEP_BUILD_DIR="BUILD"
RTOOLS_DIR=$(cygpath $FOLDER_MODE $RTOOLS_DIR)

OPENMP_LIB_FLAG=""
if [[ ! -z $SHLIB_OPENMP_CXXFLAGS ]]; then
    OPENMP_LIB_FLAG="-lgomp.dll -lgomp"
    # SHLIB_OPENMP_CXXFLAGS="-D_OPENMP $SHLIB_OPENMP_CXXFLAGS"
fi

echo "PRE CONFIGURE INJECT..."
echo "BUILD DIR: $NCBI_DIR"
echo "SCRIPTPATH: $SCRIPTPATH"
echo "RTOOLS DIR: $RTOOLS_DIR"
echo "OpenMP Flags : $SHLIB_OPENMP_CXXFLAGS"

#cd $NCBI_DIR/ncbi-cxx-toolkit-public
cd $NCBI_DIR
# export CMAKE_CMD="$RTOOLS_DIR/x86_64-w64-mingw32.static.posix/bin/cmake.exe"

# curl -OJL https://cran.r-project.org/src/base/R-latest.tar.gz
# tar xzvf R-latest.tar.gz
# cd $(find . -maxdepth 1 -iname "R-*" -type d )

# export PATH="/x86_64-w64-mingw32.static.posix/bin/:$RTOOLS_DIR/x86_64-w64-mingw32.static.posix/bin/:$PATH"
# export CC="$RTOOLS_DIR/x86_64-w64-mingw32.static.posix/bin/gcc.exe"
# export CXX="$RTOOLS_DIR/x86_64-w64-mingw32.static.posix/bin/g++.exe"
# export CMAKE_C_COMPILER=$CC
# export CMAKE_CXX_COMPILER=$CXX

#echo "GCC VER : $COMPILER_VER"
# pacman -Sy -q --needed --noconfirm development mingw-w64-cross-toolchain libraries compression mingw-w64-cross sys-utils mingw-w64-clang-x86_64-toolchain mingw-w64-ucrt-x86_64-toolchain mingw-w64-x86_64-toolchain
# pacman -Sy -q --needed --noconfirm cmake mingw-w64-clang-x86_64-cmake mingw-w64-ucrt-x86_64-cmake git mingw-w64-ucrt-x86_64-libutf8proc mingw-w64-ucrt-x86_64-utf8cpp mingw-w64-clang-x86_64-utf8cpp mingw-w64-clang-x86_64-libutf8proc mingw-w64-x86_64-libutf8proc mingw-w64-x86_64-utf8cpp mingw-w64-ucrt-x86_64-xsimd mingw-w64-x86_64-xsimd mingw-w64-clang-x86_64-xsimd mingw-w64-ucrt-x86_64-re2 mingw-w64-clang-x86_64-re2 mingw-w64-x86_64-re2 mingw-w64-ucrt-x86_64-boost mingw-w64-clang-x86_64-boost mingw-w64-x86_64-boost mingw-w64-ucrt-x86_64-thrift mingw-w64-clang-x86_64-thrift inetutils mingw-w64-ucrt-x86_64-python-win_inet_pton mingw-w64-clang-x86_64-python-win_inet_pton mingw-w64-x86_64-python-win_inet_pton mingw-w64-clang-x86_64-clang mingw-w64-ucrt-x86_64-clang clang llvm mingw-w64-clang-x86_64-llvm mingw-w64-ucrt-x86_64-llvm mingw-w64-x86_64-llvm mingw-w64-ucrt-x86_64-llvm-libs mingw-w64-clang-x86_64-llvm-libs mingw-w64-x86_64-llvm-libs lzop mingw-w64-clang-x86_64-lzo2 mingw-w64-ucrt-x86_64-lzo2 mingw-w64-x86_64-lzo2 liblzo2 liblzo2-devel mingw-w64-x86_64-openmp mingw-w64-ucrt-x86_64-openmp mingw-w64-clang-x86_64-openmp brotli mingw-w64-clang-x86_64-brotli mingw-w64-ucrt-x86_64-brotli mingw-w64-x86_64-brotli brotli-devel mingw-w64-ucrt-x86_64-samtools mingw-w64-clang-x86_64-samtools mingw-w64-x86_64-samtools mingw-w64-ucrt-x86_64-arrow mingw-w64-x86_64-arrow liblz4 lz4 mingw-w64-clang-x86_64-lz4 mingw-w64-ucrt-x86_64-lz4 mingw-w64-x86_64-lz4 zlib zlib-devel mingw-w64-cross-zlib mingw-w64-x86_64-zlib mingw-w64-ucrt-x86_64-zlib mingw-w64-clang-x86_64-zlib mingw-w64-clang-x86_64-python-conan mingw-w64-ucrt-x86_64-python-conan mingw-w64-x86_64-libc++ mingw-w64-ucrt-x86_64-libc++ mingw-w64-clang-x86_64-libc++ mingw-w64-clang-x86_64-headers-git mingw-w64-ucrt-x86_64-headers-git mingw-w64-x86_64-headers-git mingw-w64-ucrt-x86_64-libunwind mingw-w64-x86_64-libunwind mingw-w64-clang-x86_64-libunwind mingw-w64-x86_64-libbacktrace mingw-w64-ucrt-x86_64-libbacktrace mingw-w64-clang-x86_64-libbacktrace #msys2-runtime-3.5-devel msys2-runtime-3.5  msys2-w32api-runtime
# git clone https://github.com/apache/arrow
# mkdir -p $NCBI_DIR/arrow/cpp/$DEP_BUILD_DIR
# cd $NCBI_DIR/arrow/cpp/$DEP_BUILD_DIR
##export CMAKE_ARGS="-G 'MinGW Makefiles' -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS:bool=OFF -DCMAKE_POSITION_INDEPENDENT_CODE:bool=ON"
# cmake .. -G "Unix Makefiles" -DARROW_CSV=ON -DARROW_FILESYSTEM=ON -DARROW_PARQUET=ON -DARROW_THRIFT_USE_SHARED=OFF -DARROW_WITH_BZ2=ON -DARROW_BZ2_USE_SHARED=OFF -DARROW_WITH_LZ4=ON -DARROW_LZ4_USE_SHARED=OFF -DARROW_DEPENDENCY_SOURCE=AUTO -DARROW_BUILD_STATIC=ON -DARROW_BUILD_SHARED=OFF -DARROW_BOOST_USE_SHARED=OFF -DPARQUET_MINIMAL_DEPENDENCY=OFF -DARROW_BROTLI_USE_SHARED=OFF -DARROW_PROTOBUF_USE_SHARED=OFF -DARROW_SNAPPY_USE_SHARED=OFF -DARROW_UTF8PROC_USE_SHARED=OFF -DARROW_WITH_BROTLI=ON -DARROW_WITH_ZLIB=ON -DARROW_WITH_ZSTD=ON -DARROW_ZSTD_USE_SHARED=OFF -DCMAKE_BUILD_TYPE=Release -DARROW_IPC=ON -DBUILD_SHARED_LIBS:bool=OFF -DCMAKE_POSITION_INDEPENDENT_CODE:bool=ON # -DARROW_COMPUTE=ON
# make -j $(nproc)
# cd $NCBI_DIR

# git clone https://github.com/ncbi/ncbi-cxx-toolkit-public
# cp -r $(dirname $0)/src ncbi-cxx-toolkit-public/
# cp -r $(dirname $0)/include ncbi-cxx-toolkit-public/
# cp -f $(dirname $0)/QuickBLAST.prj ncbi-cxx-toolkit-public/
# mkdir -p $NCBI_DIR/ncbi-cxx-toolkit-public/$DEP_BUILD_DIR
# cd $NCBI_DIR/ncbi-cxx-toolkit-public

#export CMAKE_ARGS="-G 'MSYS Makefiles' -DCMAKE_SYSTEM_NAME=Linux -DCMAKE_SYSTEM_PROCESSOR=x86_64 -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS:bool=OFF -DCMAKE_POSITION_INDEPENDENT_CODE:bool=ON" #-DUNIX=1 -DWIN32=0

# echo "SHELL TYPE : $SHELL_TYPE"

#############GCC - OLD
# rm -f $NCBI_DIR/ncbi-cxx-toolkit-public/src/build-system/cmake/toolchains/$toolchain_file
# cp -f $NCBI_DIR/ncbi-cxx-toolkit-public/src/build-system/cmake/toolchains/x86_64-linux-gcc.cmake.in $NCBI_DIR/ncbi-cxx-toolkit-public/src/build-system/cmake/toolchains/$toolchain_file
# #sed -i "s/-Wl,--enable-new-dtags//g" $NCBI_DIR/ncbi-cxx-toolkit-public/src/build-system/cmake/toolchains/$toolchain_file #mingw gcc does not support --enable-new-dtags flag
# echo "set(CMAKE_SYSTEM_PROCESSOR x86_64)" >> $NCBI_DIR/ncbi-cxx-toolkit-public/src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_SYSTEM_NAME Linux)" >> $NCBI_DIR/ncbi-cxx-toolkit-public/src/build-system/cmake/toolchains/$toolchain_file
# echo "set(NCBI_COMPILER GNU)" >> $NCBI_DIR/ncbi-cxx-toolkit-public/src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_SYSTEM_VERSION 22.0.4)" >> $NCBI_DIR/ncbi-cxx-toolkit-public/src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_SYSROOT /$SHELL_TYPE/ /x86_64-w64-mingw32.static.posix /usr /opt)" >> $NCBI_DIR/ncbi-cxx-toolkit-public/src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_SYSROOT_COMPILE /$SHELL_TYPE/include /x86_64-w64-mingw32.static.posix/include /usr/include /usr/lib/gcc/x86_64-pc-msys/$COMPILER_VER/include /$SHELL_TYPE/include/ /usr/lib/gcc/x86_64-pc-msys/$COMPILER_VER/include /$SHELL_TYPE/include/c++  /usr/lib/gcc/x86_64-pc-msys/$COMPILER_VER/include/c++)" >> $NCBI_DIR/ncbi-cxx-toolkit-public/src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_SYSROOT_LINK /$SHELL_TYPE/lib  /x86_64-w64-mingw32.static.posix/lib /usr/lib)" >> $NCBI_DIR/ncbi-cxx-toolkit-public/src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_STAGING_PREFIX $NCBI_DIR/ncbi-cxx-toolkit-public/$DEP_BUILD_DIR/build)" >> $NCBI_DIR/ncbi-cxx-toolkit-public/src/build-system/cmake/toolchains/$toolchain_file
# # echo "set(CMAKE_C_COMPILER $CC)" >> $NCBI_DIR/ncbi-cxx-toolkit-public/src/build-system/cmake/toolchains/$toolchain_file
# # echo "set(CMAKE_CXX_COMPILER $CXX)" >> $NCBI_DIR/ncbi-cxx-toolkit-public/src/build-system/cmake/toolchains/$toolchain_file
# # echo "set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)" >> $NCBI_DIR/ncbi-cxx-toolkit-public/src/build-system/cmake/toolchains/$toolchain_file
# # echo "set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)" >> $NCBI_DIR/ncbi-cxx-toolkit-public/src/build-system/cmake/toolchains/$toolchain_file
# # echo "set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)" >> $NCBI_DIR/ncbi-cxx-toolkit-public/src/build-system/cmake/toolchains/$toolchain_file
# # echo "set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE ONLY)" >> $NCBI_DIR/ncbi-cxx-toolkit-public/src/build-system/cmake/toolchains/$toolchain_file
# echo "set(MINGW32 1)" >> $NCBI_DIR/ncbi-cxx-toolkit-public/src/build-system/cmake/toolchains/$toolchain_file
#######################################################################

###################################GCC-NEW
## export CC="gcc" #cygwin-msys - yes #linux - no #linux-msys
## export CXX="g++" #cygwin-msys - yes #linux - no #linux-msys

export CC="/usr/bin/gcc"  #linux - no
export CXX="/usr/bin/g++" #linux - no
# export CC=$(cygpath $FOLDER_MODE "/mingw64/bin/clang")
# export CXX=$(cygpath $FOLDER_MODE "/mingw64/bin/clang++")
# export CC=$(cygpath $FOLDER_MODE "/clang64/bin/clang")
# export CXX=$(cygpath $FOLDER_MODE "/clang64/bin/clang++")

#export MAKE_PRG="make" #linux - no
# export CC="$RTOOLS_DIR/mingw64/bin/gcc.exe"  #linux - no
# export CXX="$RTOOLS_DIR/mingw64/bin/g++.exe" #linux - no
# #export MAKE_PRG="/mingw64/bin/mingw32-make"  #linux - no
# export CC=$(cygpath $FOLDER_MODE "/clang64/bin/clang")
# export CXX=$(cygpath $FOLDER_MODE "/clang64/bin/clang++")
# # export MAKE_PRG="mingw32-make"
## export CC="$RTOOLS_DIR/x86_64-w64-mingw32.static.posix/bin/gcc.exe"
## export CXX="$RTOOLS_DIR/x86_64-w64-mingw32.static.posix/bin/g++.exe"

# GCC_INC_PATH=$(cygpath $FOLDER_MODE $(x86_64-pc-msys-gcc -print-search-dirs | grep -i "install" | awk -F':' '{print $2}'))
#$(x86_64-pc-msys-gcc -print-search-dirs | grep -i "install" | awk -F':' '{print $2}' | tr '[:space:]' '\0') #;

COMPILER_VER=$($CC -dumpversion | tr -d '[:space:]')
COMPILER_VER_NOPUNCT=$($CC -dumpversion | tr -d '[:space:]' | tr -d '[:punct:]')

# COMPILER_VER=$($CC --version | sed -n 's/clang version \([0-9.]\+\).*/\1/p')
# COMPILER_VER_NOPUNCT=$($CXX --version | sed -n 's/clang version \([0-9.]\+\).*/\1/p' | tr -d '[:punct:]')

toolchain_file=$(src/build-system/cmake/toolchains/cmkTool.sh GCC $COMPILER_VER | grep "in" | awk -F' ' '{print $NF}')
# toolchain_file=$(src/build-system/cmake/toolchains/cmkTool.sh Clang $COMPILER_VER | grep "in" | awk -F' ' '{print $NF}')

#toolchain_file=$1
if [[ -z $toolchain_file ]]; then
    #toolchain_file="toolchain.in"
    echo "EXITING. COULD NOT GUESS TOOLCHAIN"
    exit 0
fi

# DATATOOL_PATH=$(cygpath $FOLDER_MODE $NCBI_DIR/$DEP_BUILD_DIR/bin/)
SHELL_TYPE=$(echo $MSYSTEM_CARCH | tr "[:upper:]" "[:lower:]")
echo "TOOLCHAIN FILE : $toolchain_file"
echo "CC : $CC"
echo "CXX : $CXX"
# echo "GCC INC PATH : $GCC_INC_PATH"
echo "COMPILER VER : $COMPILER_VER"
# echo "SHELL TYPE : $SHELL_TYPE"

##########CLANG
# ln -sf /clang64/bin/clang /clang64/bin/clang-$(/clang64/bin/clang -dumpversion)
# ln -sf /clang64/bin/clang++ /clang64/bin/clang++-$(/clang64/bin/clang++ -dumpversion)
ln -sf $CC $(realpath $CC)-$COMPILER_VER
ln -sf $CXX $(realpath $CXX)-$COMPILER_VER

# ln -sf /mingw64/bin/clang-cl /mingw64/bin/clang-"$COMPILER_VER"
# ln -sf /mingw64/bin/clang-cl /mingw64/bin/clang++-"$COMPILER_VER"
##############

rm -f src/build-system/cmake/toolchains/$toolchain_file

#############GCC
# cp -f src/build-system/cmake/toolchains/x86_64-linux-gcc.cmake.in src/build-system/cmake/toolchains/$toolchain_file
###############
##############CLANG
touch src/build-system/cmake/toolchains/$toolchain_file
################

#############GCC
echo "set(NCBI_PTBCFG_FLAGS_DEFINED YES)" >>src/build-system/cmake/toolchains/$toolchain_file
echo "include_guard(GLOBAL)" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(CMAKE_C_FLAGS_INIT         \"-gdwarf-4 -Wall -Wno-format-y2k\")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(CMAKE_CXX_FLAGS_INIT       \"-gdwarf-4 -Wall -Wno-format-y2k\")" >>src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_C_FLAGS_INIT         \"-gdwarf-2 -mfpmath=sse -msse2 -mstackrealign -Wall -Wno-format-y2k\")" >>src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_CXX_FLAGS_INIT       \"-gdwarf-2 -mfpmath=sse -msse2 -mstackrealign -Wall -Wno-format-y2k\")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(NCBI_COMPILER_FLAGS_SSE       \"-msse4.2\")" >>src/build-system/cmake/toolchains/$toolchain_file
# echo "set(NCBI_COMPILER_FLAGS_SSE       \"-msse2\")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(NCBI_COMPILER_FLAGS_COVERAGE  \"--coverage\")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(NCBI_LINKER_FLAGS_COVERAGE     \"--coverage\")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(NCBI_COMPILER_FLAGS_MAXDEBUG  \"-fsanitize=address -fstack-check\")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(NCBI_LINKER_FLAGS_MAXDEBUG   \"-fsanitize=address\")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(NCBI_LINKER_FLAGS_STATICCOMPONENTS \"-static-libgcc -static-libstdc++\")" >>src/build-system/cmake/toolchains/$toolchain_file
################

##sed -i "s/-Wl,--enable-new-dtags//g" src/build-system/cmake/toolchains/$toolchain_file #mingw gcc does not support --enable-new-dtags flag
##--host=x86_64-unknown-cygwin --build=x86_64-w64-mingw32.static.posix
# echo "set(CMAKE_HOST_SYSTEM_PROCESSOR x86_64)" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(CMAKE_HOST_SYSTEM_NAME cygwin)" >>src/build-system/cmake/toolchains/$toolchain_file
# # # echo "set(CMAKE_HOST_SYSTEM_NAME linux)" >>src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_SYSTEM_PROCESSOR x86_64)" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(CMAKE_SYSTEM_NAME msys)" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(CMAKE_BUILD_TYPE Release)" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(BUILD_SHARED_LIBS ON)" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(CMAKE_POSITION_INDEPENDENT_CODE ON)" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(NCBI_COMPILER GNU)" >>src/build-system/cmake/toolchains/$toolchain_file
# echo "set(NCBI_COMPILER Clang)" >>src/build-system/cmake/toolchains/$toolchain_file

##CLANG
echo "set(NCBI_COMPILER_VERSION $COMPILER_VER_NOPUNCT)" >>src/build-system/cmake/toolchains/$toolchain_file
###

echo "set(CMAKE_STATIC_LIBRARY_SUFFIX_C \".a\")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(CMAKE_STATIC_LIBRARY_PREFIX_C \"lib\")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(CMAKE_STATIC_LIBRARY_SUFFIX_CXX \".a\")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(CMAKE_STATIC_LIBRARY_PREFIX_CXX \"lib\")" >>src/build-system/cmake/toolchains/$toolchain_file

# echo "set(CMAKE_SHARED_LIBRARY_SUFFIX_C \".dll.a\")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(CMAKE_SHARED_LIBRARY_SUFFIX_C \".dll\")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(CMAKE_SHARED_LIBRARY_PREFIX_C \"lib\")" >>src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_SHARED_LIBRARY_SUFFIX_CXX \".dll.a\")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(CMAKE_SHARED_LIBRARY_SUFFIX_CXX \".dll\")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(CMAKE_SHARED_LIBRARY_PREFIX_CXX \"lib\")" >>src/build-system/cmake/toolchains/$toolchain_file

# echo "option(BUILD_SHARED_LIBS \"Build using shared libraries\" OFF)" >> src/build-system/cmake/toolchains/$toolchain_file

# echo "set(CMAKE_SYSTEM_VERSION 22.0.4)" >> src/build-system/cmake/toolchains/$toolchain_file
#echo "set(CMAKE_SYSROOT $RTOOLS_DIR/usr $RTOOLS_DIR/opt)" >> src/build-system/cmake/toolchains/$toolchain_file #$RTOOLS_DIR/x86_64-w64-mingw32.static.posix $RTOOLS_DIR/mingw64/ $RTOOLS_DIR/ucrt64/
# echo "set(CMAKE_SYSROOT_COMPILE \"$RTOOLS_DIR/usr/include;$GCC_INC_PATH/include;\")" >>src/build-system/cmake/toolchains/$toolchain_file # $RTOOLS_DIR/x86_64-w64-mingw32.static.posix/include $RTOOLS_DIR/mingw64/include $RTOOLS_DIR/ucrt64/include

# # #echo "set(CMAKE_SYSROOT \"$RTOOLS_DIR\")" >> src/build-system/cmake/toolchains/$toolchain_file
# # echo "list(APPEND CMAKE_SYSROOT \"$RTOOLS_DIR/usr/include\")" >> src/build-system/cmake/toolchains/$toolchain_file
# # echo "list(APPEND CMAKE_SYSROOT \"$GCC_INC_PATH/include\")" >> src/build-system/cmake/toolchains/$toolchain_file
# # echo "list(APPEND CMAKE_SYSROOT \"$RTOOLS_DIR/opt\")" >> src/build-system/cmake/toolchains/$toolchain_file
# # echo "list(APPEND CMAKE_PREFIX_PATH \"$RTOOLS_DIR/usr/include\")" >> src/build-system/cmake/toolchains/$toolchain_file
# # echo "list(APPEND CMAKE_PREFIX_PATH \"$GCC_INC_PATH/include\")" >> src/build-system/cmake/toolchains/$toolchain_file
# # echo "set(CMAKE_SYSROOT_COMPILE \"$RTOOLS_DIR\")" >>src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_SYSROOT_COMPILE \"$GCC_INC_PATH/include\")" >>src/build-system/cmake/toolchains/$toolchain_file

# echo "set(CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES \"$RTOOLS_DIR/usr/\")" >>src/build-system/cmake/toolchains/$toolchain_file
# echo "list(APPEND CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES \"$GCC_INC_PATH/\")" >>src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES \"$RTOOLS_DIR/usr/\")" >>src/build-system/cmake/toolchains/$toolchain_file
# echo "list(APPEND CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES \"$GCC_INC_PATH/\")" >>src/build-system/cmake/toolchains/$toolchain_file

# echo "set(CMAKE_SYSROOT_LINK \"/x86_64-w64-mingw32.static.posix/\")" >>src/build-system/cmake/toolchains/$toolchain_file #GCC
# echo "set(CMAKE_SYSROOT_COMPILE \"$RTOOLS_DIR/mingw64/\")" >>src/build-system/cmake/toolchains/$toolchain_file           #GCC

# echo "set(CMAKE_STAGING_PREFIX $NCBI_DIR/ncbi-cxx-toolkit-public/$DEP_BUILD_DIR/build)" >> src/build-system/cmake/toolchains/$toolchain_file

#-DHAVE_IOSTREAM=1 -DNCBI_OS_OSF1=1 -DHAVE_INTTYPES_H=1 -DHAVE_NETINET_TCP_H=1 -D'_DEBUG_ARG(arg)=arg'
#-include $NCBI_DIR/ncbi_preprocessor_hooks.h
# -U_DEBUG_ARG
#-static-libgcc -static-libstdc++

echo "set(CMAKE_C_COMPILER $CC)" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(CMAKE_CXX_COMPILER $CXX)" >>src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_MAKE_PROGRAM  $MAKE_PRG)" >>src/build-system/cmake/toolchains/$toolchain_file #mingw uses different make

# #############COMPILER - MINGW64 - GCC################################# - WONT COMPILE
# mkdir -p $NCBI_DIR/$DEP_BUILD_DIR/inc/sys
# cp -f $RTOOLS_DIR/mingw64/include/errno.h $NCBI_DIR/$DEP_BUILD_DIR/inc/sys
# echo "set(CMAKE_C_FLAGS_RELEASE \"-I$RTOOLS_DIR/mingw64/include/ -idirafter $RTOOLS_DIR/mingw64/include/ -std=gnu11 -g -O2 -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_DEBUG -DHAVE_IOSTREAM=1 -DNCBI_OS_OSF1=1 -DHAVE_INTTYPES_H=1 -DHAVE_NETINET_TCP_H=1 -D_FORTIFY_SOURCE=2 -UNDEBUG -Wall -pedantic -g -O0 -fdiagnostics-color=always -fPIC -O3 $SHLIB_OPENMP_CXXFLAGS -pthread -fpermissive \")" >>src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_CXX_FLAGS_RELEASE \"-I$RTOOLS_DIR/mingw64/include/ -idirafter $RTOOLS_DIR/mingw64/include/ -std=gnu++17 -g -O2 -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_DEBUG -DHAVE_IOSTREAM=1 -DNCBI_OS_OSF1=1 -DHAVE_INTTYPES_H=1 -DHAVE_NETINET_TCP_H=1 -D_FORTIFY_SOURCE=2 -UNDEBUG -Wall -pedantic -g -O0 -fdiagnostics-color=always -fPIC -O3 $SHLIB_OPENMP_CXXFLAGS -pthread -fpermissive \")" >>src/build-system/cmake/toolchains/$toolchain_file
# ####################################################

#############COMPILER - MSYS2 - GCC#################################
#-flto #-mabi=ms -mdll
#-static-libstdc++
#-static-libgcc
# echo "set(CMAKE_C_FLAGS_RELEASE \" -v -static -std=gnu11 -O2 -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -U_WIN32 -D_DEBUG -DHAVE_IOSTREAM=1 -DNCBI_OS_OSF1=1 -DHAVE_INTTYPES_H=1 -DHAVE_NETINET_TCP_H=1 -D_FORTIFY_SOURCE=2 -UNDEBUG -Wall -pedantic -O0 -fdiagnostics-color=always -fPIC -O3 -pthread -mthreads -fvisibility=default -mtune=generic -maccumulate-outgoing-args -Wl,--allow-shlib-undefined -Wl,--unresolved-symbols=ignore-all \")" >>src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_CXX_FLAGS_RELEASE \" -v -static -std=gnu++17 -O2 -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -U_WIN32 -D_DEBUG -DHAVE_IOSTREAM=1 -DNCBI_OS_OSF1=1 -DHAVE_INTTYPES_H=1 -DHAVE_NETINET_TCP_H=1 -D_FORTIFY_SOURCE=2 -UNDEBUG -Wall -pedantic -O0 -fdiagnostics-color=always -fPIC -O3 $SHLIB_OPENMP_CXXFLAGS -pthread -mthreads -fpermissive -fvisibility=default -mtune=generic -maccumulate-outgoing-args -Wl,--allow-shlib-undefined -Wl,--unresolved-symbols=ignore-all \")" >>src/build-system/cmake/toolchains/$toolchain_file
# -Wl,--allow-shlib-undefined -Wl,--unresolved-symbols=ignore-all
#-Wl,-Bdynamic -Wl,--export-all-symbols -Wl,--demangle
#-gdwarf-2 -mfpmath=sse -msse2 -mstackrealign
#-O3
echo "set(CMAKE_C_FLAGS_RELEASE \" -std=gnu11 -O2 -Wformat -Werror=format-security -Wdate-time -U_WIN32 -D_DEBUG -DHAVE_IOSTREAM=1 -DNCBI_OS_OSF1=1 -DHAVE_INTTYPES_H=1 -DHAVE_NETINET_TCP_H=1 -D_FORTIFY_SOURCE=2 -UNDEBUG -Wall -pedantic -fstack-protector-strong -fdiagnostics-color=always -pthread -mthreads -mtune=generic -maccumulate-outgoing-args -Wall -Wno-format-y2k \")" >>src/build-system/cmake/toolchains/$toolchain_file
echo "set(CMAKE_CXX_FLAGS_RELEASE \" -std=gnu++17 -O2 -Wformat -Werror=format-security -Wdate-time -U_WIN32 -D_DEBUG -DHAVE_IOSTREAM=1 -DNCBI_OS_OSF1=1 -DHAVE_INTTYPES_H=1 -DHAVE_NETINET_TCP_H=1 -D_FORTIFY_SOURCE=2 -UNDEBUG -Wall -pedantic -fstack-protector-strong -fdiagnostics-color=always $SHLIB_OPENMP_CXXFLAGS -pthread -mthreads -mtune=generic -maccumulate-outgoing-args -Wall -Wno-format-y2k \")" >>src/build-system/cmake/toolchains/$toolchain_file
###################################################

# #############COMPILER - CLANG#################################
# # -DNCBI_OS_OSF1=1 -DHAVE_INTTYPES_H=1 -DHAVE_NETINET_TCP_H=1 -I$(cygpath $FOLDER_MODE /clang64/include/) -I$(cygpath $FOLDER_MODE /clang64/include/c++/v1) --target=x86_64-w64-windows-gnu
# # echo "set(NCBI_WIN32_THREADS 1)" >>src/build-system/cmake/toolchains/$toolchain_file
# # echo "set(NCBI_OS_MSWIN 1)" >>src/build-system/cmake/toolchains/$toolchain_file
# # echo "set(NCBI_OS_UNIX 0)" >>src/build-system/cmake/toolchains/$toolchain_file
# #--target=x86_64-unknown-windows-gnu -D_MSC_VER=$COMPILER_VER_NOPUNCT -D_WIN32 -I$(cygpath $FOLDER_MODE /clang64/include/c++/v1)
# echo "set(CMAKE_C_FLAGS \" -v -std=gnu11 -O2 -Wformat -Werror=format-security -Wdate-time -D_DEBUG -DHAVE_LIMITS=1 -DHAVE_IOSTREAM=1 -D_FORTIFY_SOURCE=2 -DWIN32 -D_WIN32 -D_WINDOWS -D_MSC_VER=$COMPILER_VER_NOPUNCT -UNDEBUG -Wall -pedantic -O0 -fdiagnostics-color=always -O3 -pthread -mtune=generic  \")" >>src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_CXX_FLAGS \" -v -std=gnu++17 -O2 -Wformat -Werror=format-security -Wdate-time -D_DEBUG -DHAVE_LIMITS=1 -DHAVE_IOSTREAM=1 -D_FORTIFY_SOURCE=2 -DWIN32 -D_WIN32 -D_WINDOWS -D_MSC_VER=$COMPILER_VER_NOPUNCT -UNDEBUG -Wall -pedantic -O0 -fdiagnostics-color=always -O3 -D_OPENMP -fopenmp -pthread -mtune=generic \")" >>src/build-system/cmake/toolchains/$toolchain_file
# # ####################################################

##############LINKER - GCC##############
#-Wl,--architecture=x86_64 -m64 -L/usr/lib/ $OPENMP_LIB_FLAG -Wl,-ffat-lto-objects #-Wl,--allow-shlib-undefined -Wl,--unresolved-symbols=ignore-all -Wl,--warn-unresolved-symbols #-Wl,--script=$SCRIPTPATH/arrow_int64_t_linker.ld #-Wl,--no-undefined
#-Wl,--compat-implib
echo "set(CMAKE_EXE_LINKER_FLAGS_INIT  \" -Wl,--demangle -Wl,--pic-executable -Wl,--support-old-code -Wl,--as-needed \")" >>src/build-system/cmake/toolchains/$toolchain_file
#-Wl,-ffat-lto-objects -Wl,--allow-shlib-undefined -Wl,--unresolved-symbols=ignore-all -Wl,-flto=auto  -Wl,--dll -Wl,-fstack-protector-strong
echo "set(CMAKE_SHARED_LINKER_FLAGS_INIT  \" -Wl,--demangle -Wl,--pic-executable -Wl,-ffat-lto-objects -Wl,-flto=auto -Wl,-fstack-protector-strong -Wl,--enable-auto-import -Wl,--support-old-code -Wl,--as-needed -Wl,--dll -Wl,-shared \")" >>src/build-system/cmake/toolchains/$toolchain_file
################################

# ###################LINKER - CLANG#####################
# # echo "set(CMAKE_EXE_LINKER_FLAGS_INIT  \"-L/clang64/lib/ -L/x86_64-w64-mingw32.static.posix/lib/ -Wl,--demangle -Wl,--pic-executable -Wl,--export-all-symbols -Wl,--as-needed\")" >>src/build-system/cmake/toolchains/$toolchain_file
# # echo "set(CMAKE_SHARED_LINKER_FLAGS_INIT  \"-L/clang64/lib/ -L/x86_64-w64-mingw32.static.posix/lib/ -Wl,--demangle -Wl,--pic-executable -Wl,-shared -Wl,--export-all-symbols -Wl,-flto=auto -Wl,-ffat-lto-objects -Wl,--no-undefined -Wl,--as-needed\")" >>src/build-system/cmake/toolchains/$toolchain_file
# # -Wl,-shared -Wl,--dll
# echo "set(CMAKE_EXE_LINKER_FLAGS_INIT  \" -Wl,-v -Wl,--demangle -Wl,--pic-executable -Wl,--as-needed  \")" >>src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_SHARED_LINKER_FLAGS_INIT  \" -Wl,-v -Wl,-static -Wl,--demangle -Wl,--pic-executable -Wl,--as-needed \")" >>src/build-system/cmake/toolchains/$toolchain_file
# # #####################################

# echo "set(CMAKE_LIBRARY_PATH \"/x86_64-w64-mingw32.static.posix/lib/;/usr/lib/;mingw64/lib/\")" >>src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_PROGRAM_PATH \"/x86_64-w64-mingw32.static.posix/bin/;/usr/bin/;/mingw64/bin/\")" >>src/build-system/cmake/toolchains/$toolchain_file

echo "include($NCBI_DIR/ncbi_components.cmake)" >>src/build-system/cmake/toolchains/$toolchain_file

# echo "message(\"NCBI INC ROOT : \" $NCBI_DIR/include)" >>src/build-system/cmake/toolchains/$toolchain_file

# echo "set(PATH \$ENV{PATH})"  >> src/build-system/cmake/toolchains/$toolchain_file
# echo "set(ENV{PATH} \"$NCBI_DIR/ncbi-cxx-toolkit-public/$DEP_BUILD_DIR/bin/:\"\${PATH})"  >> src/build-system/cmake/toolchains/$toolchain_file
# echo "message(\"PATH :\" \${PATH})"  >> src/build-system/cmake/toolchains/$toolchain_file
# echo "list(APPEND CMAKE_PROGRAM_PATH \"$NCBI_DIR/ncbi-cxx-toolkit-public/$DEP_BUILD_DIR/bin/\")" >>src/build-system/cmake/toolchains/$toolchain_file

# cat << EOF >> src/build-system/cmake/toolchains/$toolchain_file
# if(NOT DEFINED NCBI_DATATOOL)
#     set(NCBI_DATATOOL "$DATATOOL_PATH/datatool.exe")
#     FILE(STRINGS "$NCBI_DIR/src/build-system/datatool_version.txt" _datatool_version)
#     set(NCBI_DATATOOL_PATH "$DATATOOL_PATH")
#     set(NCBI_DATATOOL_BASE "$DATATOOL_PATH")
#     set(NCBI_DATATOOL_BIN "datatool.exe")
#     include($NCBI_DIR/src/build-system/cmake/CMake.NCBIptb.cmake)
#     include($NCBI_DIR/src/build-system/cmake/CMake.NCBIptb.datatool.cmake)
# endif(NOT DEFINED NCBI_DATATOOL)
# EOF

# echo "set(NCBI_DATATOOL_PATH \"$DATATOOL_PATH\")" >> src/build-system/cmake/toolchains/$toolchain_file
# echo "set(NCBI_DATATOOL_BASE \"$DATATOOL_PATH/\")" >> src/build-system/cmake/toolchains/$toolchain_file
# echo "set(NCBI_DATATOOL_BIN \"datatool.exe\")" >> src/build-system/cmake/toolchains/$toolchain_file
# echo "FILE(STRINGS \"$NCBI_DIR/src/build-system/datatool_version.txt\" _datatool_version)" >> src/build-system/cmake/toolchains/$toolchain_file
# echo "set(NCBI_DATATOOL \"$DATATOOL_PATH/datatool.exe\")" >> src/build-system/cmake/toolchains/$toolchain_file

# echo "set(NCBI_DATATOOL_PATH \"$DATATOOL_PATH\")" >>src/build-system/cmake/toolchains/$toolchain_file
# echo "set(NCBI_DATATOOL_BASE \"$DATATOOL_PATH/bin/\")" >>src/build-system/cmake/toolchains/$toolchain_file
# echo "set(NCBI_DATATOOL_BIN \"datatool.exe\")" >>src/build-system/cmake/toolchains/$toolchain_file

# echo "set(CMAKE_STATIC_LINKER_FLAGS_INIT  \"-Wl,--pic-executable -Wl,-static -Wl,-Bstatic -Wl,--support-old-code -Wl,--compat-implib -Wl,--architecture=x86_64 -Wl,--export-all-symbols -Wl,-flto=auto -Wl,-ffat-lto-objects -Wl,--no-undefined  -Wl,--as-needed\")" >> src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_EXE_LINKER_FLAGS \"--we3 --export-all-symbols -flto=auto -ffat-lto-objects \")" >> src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_SHARED_LINKER_FLAGS \"--we2 --export-all-symbols -flto=auto -ffat-lto-objects \")" >> src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_STATIC_LINKER_FLAGS \"--we1 --export-all-symbols -flto=auto -ffat-lto-objects \")" >> src/build-system/cmake/toolchains/$toolchain_file
echo "set(CMAKE_CXX_EXTENSIONS ON)" >>src/build-system/cmake/toolchains/$toolchain_file
###TRY ADDING CLANG/18/LIB/INNCLUDE

#-v -H
# -fuse-ld=lld -- is risky for msys and will not work
#-Wl,--enable-new-dtags -- is risky for msys and will not work

#-B $GCC_INC_PATH - wont work during cmake config but works during make

#-std=gnu17 -std=gnu++17
#-isystem $RTOOLS_DIR/usr/include/cygwin
# echo "set(CMAKE_INCLUDE_DIRECTORIES_BEFORE ON)" >> src/build-system/cmake/toolchains/$toolchain_file
#-I $RTOOLS_DIR/usr/include -I $GCC_INC_PATH/include/c++ -I $GCC_INC_PATH/include/c++/x86_64-pc-msys/
#-isysroot $RTOOLS_DIR/clang64/include/ -I $RTOOLS_DIR/clang64/lib/clang/18/include -isystem $RTOOLS_DIR/clang64/include/c++/v1 -idirafter $RTOOLS_DIR/clang64/include

#############TESTING
# echo "set(CMAKE_C_FLAGS_RELEASE \"-v -H -nostdinc -g -O2 -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -DHAVE_IOSTREAM -D_FORTIFY_SOURCE=2  -UNDEBUG -Wall -pedantic -g -O0 -fdiagnostics-color=always -fPIC -O3 -fuse-ld=lld -fpermissive -isystem $RTOOLS_DIR/usr/include/ -isystem $GCC_INC_PATH/include/ \")" >> src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_CXX_FLAGS_RELEASE \"-v -H -nostdinc -g -O2 -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -DHAVE_IOSTREAM -D_FORTIFY_SOURCE=2  -UNDEBUG -Wall -pedantic -g -O0 -fdiagnostics-color=always -fPIC -O3 $SHLIB_OPENMP_CXXFLAGS -pthread -fuse-ld=lld -fpermissive -isystem $GCC_INC_PATH -isystem $GCC_INC_PATH/include/c++/ -idirafter $GCC_INC_PATH/include/c++/x86_64-pc-msys/ -idirafter $RTOOLS_DIR/usr/include -isystem $GCC_INC_PATH/include/ \")" >> src/build-system/cmake/toolchains/$toolchain_file
########################

#-isystem $GCC_INC_PATH/include/c++/ -idirafter $GCC_INC_PATH/include/c++/x86_64-pc-msys/
#-I $RTOOLS_DIR/usr/include -isystem $GCC_INC_PATH -I $GCC_INC_PATH/include/c++/x86_64-pc-msys/ -isystem $GCC_INC_PATH/include/c++/x86_64-pc-msys/ -I $GCC_INC_PATH/include/c++/ -isystem $GCC_INC_PATH/include/c++/ -I $GCC_INC_PATH/include/ -isystem $GCC_INC_PATH/include/
#-nostdinc
#-Wl,-Bsymbolic-functions -Wl,-z,relro
#-isystem $GCC_INC_PATH/include/ -isystem $GCC_INC_PATH/include/c++/ -isystem $GCC_INC_PATH/include/c++/x86_64-pc-msys
# echo "set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)" >> $NCBI_DIR/ncbi-cxx-toolkit-public/src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)" >> $NCBI_DIR/ncbi-cxx-toolkit-public/src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)" >> $NCBI_DIR/ncbi-cxx-toolkit-public/src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE ONLY)" >> $NCBI_DIR/ncbi-cxx-toolkit-public/src/build-system/cmake/toolchains/$toolchain_file
echo "set(MINGW32 1)" >>src/build-system/cmake/toolchains/$toolchain_file
##################################################################################################

###############################CLANG-CL

# #Clang-cl - windows
# export CC="$RTOOLS_DIR/usr/bin/clang-cl.exe"
# export CXX="$RTOOLS_DIR/usr/bin/clang-cl.exe"
# CLANG_CL_VER=$($CC --version | head -n 1 | awk '{print $3}')
# COMPILER_VER=$(/x86_64-w64-mingw32.static.posix/bin/gcc -dumpversion | tr -d '[:space:]')

# toolchain_file=$(src/build-system/cmake/toolchains/cmkTool.sh Clang $CLANG_CL_VER | grep "in" | awk -F' ' '{print $NF}')
# toolchain_file=$1
# if [[ -z $toolchain_file ]]; then
#     #toolchain_file="toolchain.in"
#     echo "EXITING. COULD NOT GUESS TOOLCHAIN"
#     exit 0
# fi

# SHELL_TYPE=$(echo $MSYSTEM | tr "[:upper:]" "[:lower:]")
# echo "TOOLCHAIN FILE : $toolchain_file"

# WIN_VER=$(systeminfo | grep -i "OS Version" | awk '{print $3}')
# rm -f src/build-system/cmake/toolchains/$toolchain_file
# cp -f src/build-system/cmake/toolchains/amd64-windows-msvc.cmake src/build-system/cmake/toolchains/$toolchain_file
# #sed -i "s/-Wl,--enable-new-dtags//g" src/build-system/cmake/toolchains/$toolchain_file #mingw gcc does not support --enable-new-dtags flag
# echo "set(CMAKE_SYSTEM_PROCESSOR amd64)" >> src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_SYSTEM_NAME Windows)" >> src/build-system/cmake/toolchains/$toolchain_file
# echo "set(NCBI_COMPILER Clang)" >> src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_SYSTEM_VERSION $WIN_VER)" >> src/build-system/cmake/toolchains/$toolchain_file
# # echo "set(CMAKE_SYSROOT /x86_64-w64-mingw32.static.posix /usr /mingw64/ /ucrt64/ /opt)" >> src/build-system/cmake/toolchains/$toolchain_file
# # echo "set(CMAKE_SYSROOT_COMPILE /x86_64-w64-mingw32.static.posix/include /usr/include /mingw64/include /ucrt64/include)" >>src/build-system/cmake/toolchains/$toolchain_file
# # echo "set(CMAKE_SYSROOT_LINK /x86_64-w64-mingw32.static.posix/lib /mingw64/lib /ucrt64/lib /usr/lib)" >> src/build-system/cmake/toolchains/$toolchain_file
# # echo "set(CMAKE_STAGING_PREFIX $NCBI_DIR/ncbi-cxx-toolkit-public/$DEP_BUILD_DIR/build)" >> src/build-system/cmake/toolchains/$toolchain_file
# echo "message(DBG here)" >> src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_C_COMPILER $CC)" >> src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_CXX_COMPILER $CXX)" >> src/build-system/cmake/toolchains/$toolchain_file
# # echo "set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)" >> $NCBI_DIR/ncbi-cxx-toolkit-public/src/build-system/cmake/toolchains/$toolchain_file
# # echo "set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)" >> $NCBI_DIR/ncbi-cxx-toolkit-public/src/build-system/cmake/toolchains/$toolchain_file
# # echo "set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)" >> $NCBI_DIR/ncbi-cxx-toolkit-public/src/build-system/cmake/toolchains/$toolchain_file
# # echo "set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE ONLY)" >> $NCBI_DIR/ncbi-cxx-toolkit-public/src/build-system/cmake/toolchains/$toolchain_file
# echo "set(MINGW32 1)" >> src/build-system/cmake/toolchains/$toolchain_file
# ##############################################################

# ######################Clang
# export CC="/mingw64/bin/clang"
# export CXX="/mingw64/bin/clang++"
# # export CC="clang-cl"
# # export CXX="clang-cl"
# CLANG_VER=$($CC -dumpversion | tr -d '[:space:]')
# COMPILER_VER=$(/x86_64-w64-mingw32.static.posix/bin/gcc -dumpversion | tr -d '[:space:]')
# # toolchain_file=$(src/build-system/cmake/toolchains/cmkTool.sh Clang $CLANG_VER | grep "in" | awk -F' ' '{print $NF}')
# toolchain_file=$(src/build-system/cmake/toolchains/cmkTool.sh GCC $COMPILER_VER | grep "in" | awk -F' ' '{print $NF}')
# #toolchain_file=$1
# if [[ -z $toolchain_file ]]; then
#     #toolchain_file="toolchain.in"
#     echo "EXITING. COULD NOT GUESS TOOLCHAIN"
#     exit 0
# fi

# SHELL_TYPE=$(echo $MSYSTEM | tr "[:upper:]" "[:lower:]")
# echo "TOOLCHAIN FILE : $toolchain_file"
# rm -f src/build-system/cmake/toolchains/$toolchain_file
# cp -f src/build-system/cmake/toolchains/x86_64-darwin-clang.cmake.in src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_SYSTEM_PROCESSOR x86_64)" >> src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_SYSTEM_NAME cygwin)" >> src/build-system/cmake/toolchains/$toolchain_file
# echo "set(NCBI_COMPILER GNU)" >> src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_SYSTEM_VERSION 22.0.4)" >> src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_SYSROOT $RTOOLS_DIR/x86_64-w64-mingw32.static.posix $RTOOLS_DIR/usr $RTOOLS_DIR/mingw64/ $RTOOLS_DIR/ucrt64/ $RTOOLS_DIR/opt)" >> src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_SYSROOT_COMPILE $RTOOLS_DIR/x86_64-w64-mingw32.static.posix/include $RTOOLS_DIR/usr/include $RTOOLS_DIR/mingw64/include $RTOOLS_DIR/ucrt64/include)" >>src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_SYSROOT_LINK $RTOOLS_DIR/x86_64-w64-mingw32.static.posix/lib $RTOOLS_DIR/mingw64/lib $RTOOLS_DIR/ucrt64/lib $RTOOLS_DIR/usr/lib)" >> src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_STAGING_PREFIX $NCBI_DIR/ncbi-cxx-toolkit-public/$DEP_BUILD_DIR/build)" >> src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_C_COMPILER $CC)" >> src/build-system/cmake/toolchains/$toolchain_file
# echo "set(CMAKE_CXX_COMPILER $CXX)" >> src/build-system/cmake/toolchains/$toolchain_file
# # echo "set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)" >> $NCBI_DIR/ncbi-cxx-toolkit-public/src/build-system/cmake/toolchains/$toolchain_file
# # echo "set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)" >> $NCBI_DIR/ncbi-cxx-toolkit-public/src/build-system/cmake/toolchains/$toolchain_file
# # echo "set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)" >> $NCBI_DIR/ncbi-cxx-toolkit-public/src/build-system/cmake/toolchains/$toolchain_file
# # echo "set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE ONLY)" >> $NCBI_DIR/ncbi-cxx-toolkit-public/src/build-system/cmake/toolchains/$toolchain_file
# echo "set(MINGW32 1)" >> src/build-system/cmake/toolchains/$toolchain_file
# ###################################################################################################

# export PATH="/x86_64-w64-mingw32.static.posix/bin/:$PATH"

# ################CYGWin
# mkdir -p $NCBI_DIR/cygwin/

# # (cd $NCBI_DIR/cygwin/; \
# # git clone https://github.com/openSUSE/libsolv; \
# # cd libsolv; \
# # mkdir build; \
# # cd build ; \
# # cmake .. ; \
# # make; \
# # make install)

# # export LIBSOLV_LIBS=$(pkg-config  --libs /usr/local/lib/pkgconfig/libsolv.pc)
# # export LIBSOLV_CFLAGS=$(pkg-config  --cflags /usr/local/lib/pkgconfig/libsolv.pc)
# # export LIBGCRYPT_CFLAGS=$(libgcrypt-config --cflags)
# # export LIBGCRYPT_LIBS=$(libgcrypt-config --libs)
# # export SYSROOT="/usr"
# # OLD_PATH=$PATH

# #export PATH="/x86_64-w64-mingw32.static.posix/bin/:$PATH"
# curl -L --url "https://www.cygwin.com/setup-x86_64.exe"  --ssl-reqd -o $NCBI_DIR/cygwin/setup-x86_64.exe

# # (cd $NCBI_DIR/cygwin/; \
# # git clone git://cygwin.com/git/cygwin-apps/setup.git; \
# # cd setup; \
# # ./bootstrap.sh --host=x86_64-w64-mingw32 --with-libgcrypt-prefix=$(libgcrypt-config --prefix) --with-gnu-ld ; \
# # make; \
# # make strip; \
# # make distclean)
# # #--host=x86_64-pc-msys
# # #--host=x86_64-pc-msys --build=x86_64-w64-mingw32.static.posix #--with-gnu-ld

# (cd $NCBI_DIR/cygwin/  &&  powershell.exe -c ./setup-x86_64.exe -q -D -g --enable-old-keys -w  --no-write-registry -N -n -d -B -I -f -A -o  --allow-unsupported-windows -a x86_64 -s "https://mirrors.cicku.me/cygwin" -W  -C Devel,Libs )
# #-q --quiet-mode #-l . --root .
# ###############

# # export CMAKE_ARGS="-G 'Unix Makefiles' $CMAKE_ARGS -DCMAKE_TOOLCHAIN_FILE=$NCBI_DIR/ncbi-cxx-toolkit-public/src/build-system/cmake/toolchains/$toolchain_file"
# ##mv cmake-configure cmake-configure.sh
# # mkdir -p $NCBI_DIR/ncbi-cxx-toolkit-public/$DEP_BUILD_DIR/build
# # cp -f $NCBI_DIR/ncbi-cxx-toolkit-public/src/build-system/cmake/toolchains/$toolchain_file $NCBI_DIR/ncbi-cxx-toolkit-public/$DEP_BUILD_DIR/build/ncbi_toolchain.cmake
# # export C_INCLUDE_PATH="/$SHELL_TYPE/include:/$SHELL_TYPE/include/:/usr/lib/gcc/x86_64-pc-msys/$COMPILER_VER/include:/$SHELL_TYPE/include/c++:/usr/lib/gcc/x86_64-pc-msys/$COMPILER_VER/include/c++:/x86_64-w64-mingw32.static.posix/include:/usr/include"
# # export CPLUS_INCLUDE_PATH="/$SHELL_TYPE/include:/$SHELL_TYPE/include/:/usr/lib/gcc/x86_64-pc-msys/$COMPILER_VER/include:/$SHELL_TYPE/include/c++:/usr/lib/gcc/x86_64-pc-msys/$COMPILER_VER/include/c++:/x86_64-w64-mingw32.static.posix/include:/usr/include"
# echo "$CMAKE_ARGS"
# #TODO : Successfully find headers
# #MINGW32=1 CMAKE_ARGS="$CMAKE_ARGS"
# # ./cmake-configure GCC $COMPILER_VER --with-projects=QuickBLAST.prj --without-debug --with-build-root=$DEP_BUILD_DIR #--with-conan
# #--with-dll
# # export CPPFLAGS="-I/usr/include -I/mingw64/include" ;  #-I/x86_64-w64-mingw32.static.posix/include
# # export CPATH="$(CPATH)" ;
# # export C_INCLUDE_PATH="$(C_INCLUDE_PATH)" ;
# # export CPLUS_INCLUDE_PATH="$(CPLUS_INCLUDE_PATH)" ;
# # export LDFLAGS="$FOLDER_MODE x64" ;
# # which gcc
# export CC="/usr/bin/gcc" ;
# export CPP="/usr/bin/gcc -E -v -H" ;
# export CFLAGS="-v -H -g -O2 -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2  -UNDEBUG -Wall -pedantic -g -O0 -fdiagnostics-color=always -fPIC -O3 -O2" ; #-fuse-ld=/clang64/bin/lld-link #-isystem /usr/include -isystem /x86_64-w64-mingw32.static.posix/include -isystem /clang64/include -isystem /mingw64/include
# export CXX="/usr/bin/g++" ;  #-nobuiltininc -nostdinc++
# export CXXPP="/usr/bin/g++ -E -v -H" ;  #-nobuiltininc -nostdinc++
# export CXXFLAGS="-v -H -g -O2 -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2  -UNDEBUG -Wall -pedantic -g -O0 -fdiagnostics-color=always -fPIC -O3 $SHLIB_OPENMP_CXXFLAGS -static -pthread" ; #-static-openmp -stdlib=libc++ # -fuse-ld=/clang64/bin/lld-link #-cxx-isystem /usr/include -cxx-isystem /x86_64-w64-mingw32.static.posix/include -cxx-isystem /clang64/include -cxx-isystem /mingw64/include
# ./configure --host=x86_64-unknown-cygwin --build=x86_64-w64-mingw32.static.posix --without-debug --with-static --with-static-exe --with-bin-release --with-lfs --with-experimental=C++20 --without-boost-tag --with-build-root=$DEP_BUILD_DIR --with-projects=QuickBLAST.prj ; #--with-openmp
# #-DUNIX=1 -DWIN32=0 -DCMAKE_SYSTEM_NAME=Linux -DCMAKE_SYSTEM_PROCESSOR=x86_64 -DNCBI_COMPILER=gcc #-DCMAKE_SYSTEM_NAME=mingw64
# ./cmake-configure GCC $COMPILER_VER --with-projects=QuickBLAST.prj --without-debug --with-build-root=$DEP_BUILD_DIR
# cd $NCBI_DIR/ncbi-cxx-toolkit-public/$DEP_BUILD_DIR/build
# cmake .
# make -j $(nproc)

#rm -rf $NCBI_DIR/ncbi-cxx-toolkit-public
# rm -f $NCBI_DIR/ncbi-cxx-toolkit-public/src/build-system/cmake/toolchains/$toolchain_file
# read -r endline
