#!/bin/bash

#Get script directory from mingw shell for sysroot
SCRIPTPATH="$(
    cd -- "$(dirname "$0")" >/dev/null 2>&1
    pwd -P
)"
BUILD_DIR=$(cygpath -u $(pwd))
DEP_BUILD_DIR="BUILD"
RTOOLS_DIR=$(cygpath -u $RTOOLS_DIR)
TRIPLET="x86_64-pc-msys"


#find $(dirname $(which cygpath)) -type f -iname "*x86_64-pc-*"
#find $(dirname $(which cygpath)) -type f -iname "*gcc*"

if [[ $(which gcc) ]]; then 
  GCC_INC_PATH=$(cygpath -m $(realpath $(gcc -print-search-dirs | grep -i "install" | sed 's/install://g')))
elif [[ $(which "$TRIPLET"-gcc) ]]; then 
  GCC_INC_PATH=$(cygpath -m $(realpath $("$TRIPLET"-gcc -print-search-dirs | grep -i "install" | sed 's/install://g')))
elif [[ $(which x86_64-pc-cygwin-gcc) ]]; then 
  GCC_INC_PATH=$(cygpath -m $(realpath $(x86_64-pc-cygwin-gcc -print-search-dirs | grep -i "install" | sed 's/install://g')))
fi

SHELL_TYPE=$(echo $MSYSTEM | tr "[:upper:]" "[:lower:]")

echo "POST CONFIGURE INJECT..."
echo "BUILD DIR: $BUILD_DIR"
echo "SCRIPTPATH: $SCRIPTPATH"
echo "RTOOLS DIR: $RTOOLS_DIR"

cd $BUILD_DIR

#cp -f $BUILD_DIR/BUILD/inc/common/config/ncbiconf_msvc.h $BUILD_DIR/BUILD/inc/common/config/ncbiconf_unix.h