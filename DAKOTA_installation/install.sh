#!/bin/bash 

# DAKOTA path
export DAK_SRC=/Users/imaity/codes/dakota-6.14.0.src

mkdir -p $DAK_SRC/build
export DAK_BUILD=$DAK_SRC/build

cp $DAK_SRC/cmake/BuildDakotaTemplate.cmake \
        $DAK_SRC/cmake/BuildDakotaCustom.cmake

cd $DAK_BUILD

cmake -C $DAK_SRC/cmake/BuildDakotaCustom.cmake $DAK_SRC

# Add the following lines to the CMakeCache.txt file 
# in the build directory. 
#---------------------------
#CMAKE_CXX_FLAGS:STRING=-Wno-implicit-function-declaration
#CMAKE_C_FLAGS:STRING=-Wno-implicit-function-declaration
#CMAKE_Fortran_FLAGS:STRING=-fallow-argument-mismatch
#---------------------------
make 
make install

