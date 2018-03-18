#!/bin/bash

srcName=ez_pca

if [ -e $srcName.exe ]; then
    rm $srcName.exe
fi

g++ -std=c++11 $srcName.cpp -o $srcName.exe