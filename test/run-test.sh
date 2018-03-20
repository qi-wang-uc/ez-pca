#!/bin/bash

srcName=ez-pca

if [ -e $srcName ]; then
    rm $srcName
fi

cp ../bin/$srcName .

./$srcName input.txt