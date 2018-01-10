#!/bin/sh

if [ "$#" == 0 ]; then
	BUILD_TARGET="Debug"
else
	BUILD_TARGET="$1"
fi

cd build
cmake -DCMAKE_BUILD_TYPE="$BUILD_TARGET" "${@:1}" ..
make -j3
ln -sf ../res res
