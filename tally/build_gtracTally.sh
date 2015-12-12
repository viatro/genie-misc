#!/bin/bash
base_name="gtracTally"
source_file="$base_name.cc"

echo "Compiling $source_file"
g++5.1.0 -Wall -O3 -std=c++11 \
`root-config --cflags --glibs` \
-lEG \
-o $base_name \
$source_file