#!/bin/bash

#pacman -Syuu --noconfirm

#Get script directory from mingw shell for sysroot
SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
# BUILD_DIR=`pwd`
BUILD_DIR=$(cygpath -m $(pwd))
DEP_BUILD_DIR="BUILD"
RTOOLS_DIR=$(cygpath -m $RTOOLS_DIR)
# TRIPLET="x86_64-w64-mingw32.static.posix"
TRIPLET="x86_64-pc-msys"

# GCC_INC_PATH=$(cygpath -m $(x86_64-pc-msys-gcc -print-search-dirs | grep -i "install" | awk -F':' '{print $2}'))
GCC_INC_PATH=$(cygpath -m $(realpath $("$TRIPLET"-gcc -print-search-dirs | grep -i "install" |  sed 's/install://g')))
SHELL_TYPE=$(echo $MSYSTEM | tr "[:upper:]" "[:lower:]")

echo "POST CONFIGURE INJECT..."
echo "BUILD DIR: $BUILD_DIR"
echo "SCRIPTPATH: $SCRIPTPATH"
echo "RTOOLS DIR: $RTOOLS_DIR"
#cd $BUILD_DIR/ncbi-cxx-toolkit-public
cd $BUILD_DIR

# export CC="$TRIPLET"-gcc
# export CXX="$TRIPLET"-g++
# echo "set(CMAKE_C_COMPILER $CC)" >> BUILD/build/ncbi_toolchain.cmake
# echo "set(CMAKE_CXX_COMPILER $CXX)" >> BUILD/build/ncbi_toolchain.cmake

# # echo "set(CMAKE_C_FLAGS_RELEASE \" \${CMAKE_C_FLAGS_RELEASE} -nostdinc -B $GCC_INC_PATH -isysroot $RTOOLS_DIR/usr/include/ -isystem $RTOOLS_DIR/usr/include/ -isystem $GCC_INC_PATH/include/ \")" >> BUILD/build/ncbi_toolchain.cmake
# # echo "set(CMAKE_CXX_FLAGS_RELEASE \" \${CMAKE_CXX_FLAGS_RELEASE} -nostdinc -nostdinc++ -B $GCC_INC_PATH -isysroot $RTOOLS_DIR/usr/include/ -isystem $GCC_INC_PATH/include/ -isystem $GCC_INC_PATH/include/c++/x86_64-pc-msys/ -idirafter $GCC_INC_PATH/include/c++/ -idirafter $RTOOLS_DIR/usr/include/ \")" >> BUILD/build/ncbi_toolchain.cmake

# # echo "set(CMAKE_C_FLAGS_RELEASE \" \${CMAKE_C_FLAGS_RELEASE} -nostdinc -I $RTOOLS_DIR/usr/include/ -I $GCC_INC_PATH/include/ \")" >> BUILD/build/ncbi_toolchain.cmake
# # echo "set(CMAKE_CXX_FLAGS_RELEASE \" \${CMAKE_CXX_FLAGS_RELEASE} -nostdinc -nostdinc++ -I $RTOOLS_DIR/$TRIPLET/include/ -I $GCC_INC_PATH/include/ -I $GCC_INC_PATH/include/c++/ -idirafter $RTOOLS_DIR/usr/include/ -I $GCC_INC_PATH/include/c++/$TRIPLET/ \")" >> BUILD/build/ncbi_toolchain.cmake
# # sed -i "s/NCBI_OS_CYGWIN/NCBI_OS_LINUX/g" BUILD/inc/ncbiconf_unix.h
# sed -i "s/NCBI_OS_UNIX/NCBI_OS_MSWIN/g" BUILD/inc/ncbiconf_unix.h

#-B $GCC_INC_PATH - wont work during cmake config but works during make
#TODO: try mingw64 or clang too