#!/usr/bin/env bash
# All content Copyright (C) 2018 Genomics plc
set -e -u -x -o pipefail

source "$( dirname "${BASH_SOURCE[0]}" )/../activate"

"$PROJECT_SCRIPTS/get-managed-dependencies.sh"

# TODO: clean up temp directory even if errors happen

SOURCE_DIR="$PROJECT_HOME/cpp"
INSTALL_DIR="$PROJECT_HOME"
TMP_DIR="${WECALL_BUILD}/wecall-release"
mkdir -p "$TMP_DIR"
cd "$TMP_DIR"

# crap C++ compiler detection:
if [ "X" != "X$(which clang++-3.9)" ]; then
    CMAKE_CXX_COMPILER="$(which clang++-3.9)"
elif [ "X" != "X$(which clang++)" ]; then
    CMAKE_CXX_COMPILER="$(which clang++)"
else
    echo 'Failed to detect C++ compiler.'
    exit 1
fi

# crap C++ compiler detection:
if [ "X" != "X$(which gcc-4.9)" ]; then
    CMAKE_C_COMPILER="$(which gcc-4.9)"
elif [ "X" != "X$(which gcc-4.8)" ]; then
    CMAKE_C_COMPILER="$(which gcc-4.8)"
elif [ "X" != "X$(which gcc-4)" ]; then
    CMAKE_C_COMPILER="$(which gcc-4)"
elif [ "X" != "X$(which gcc)" ]; then
    CMAKE_C_COMPILER="$(which gcc)"
else
    echo 'Failed to detect C++ compiler.'
    exit 1
fi

echo "using C++ $CMAKE_CXX_COMPILER"
echo "using C $CMAKE_C_COMPILER"

# Copy all libs from ${WECALL_BUILD}/dependencies to INSTALL_DIR/lib
mkdir -p "$PROJECT_HOME/bin"
mkdir -p "$PROJECT_HOME/lib"
cp "$WECALL_BUILD"/dependencies/*/bin/* "$PROJECT_HOME/bin"
cp "$WECALL_BUILD"/dependencies/*/lib/* "$PROJECT_HOME/lib"

echo "Building..."
cmake -DCMAKE_CXX_COMPILER="$CMAKE_CXX_COMPILER" -DCMAKE_C_COMPILER="$CMAKE_C_COMPILER" -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR" -DPLATFORM="$($PROJECT_SCRIPTS/platform)" "$SOURCE_DIR"
cmake --build . --target install -- -j "$(nproc)"
