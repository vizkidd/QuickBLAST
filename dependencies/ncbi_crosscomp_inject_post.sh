#!/bin/bash

#Get script directory from mingw shell for sysroot
SCRIPTPATH="$(
    cd -- "$(dirname "$0")" >/dev/null 2>&1
    pwd -P
)"
BUILD_DIR=$(cygpath -m $(pwd))
DEP_BUILD_DIR="BUILD"
RTOOLS_DIR=$(cygpath -m $RTOOLS_DIR)
TRIPLET="x86_64-pc-msys"

GCC_INC_PATH=$(cygpath -m $(realpath $("$TRIPLET"-gcc -print-search-dirs | grep -i "install" | sed 's/install://g')))
SHELL_TYPE=$(echo $MSYSTEM | tr "[:upper:]" "[:lower:]")

echo "POST CONFIGURE INJECT..."
echo "BUILD DIR: $BUILD_DIR"
echo "SCRIPTPATH: $SCRIPTPATH"
echo "RTOOLS DIR: $RTOOLS_DIR"

cd $BUILD_DIR
