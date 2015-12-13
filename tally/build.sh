#!/bin/bash

[[ $# == 1 ]] && ( [ "$1" == "egst" ] || [ "$1" == "gst" ] || [ "$1" == "gtrac" ] ) || ( echo "$0 [egst|gst|gtrac]" && exit )

base_name="$1"
source_file="$1Tally.cc"
host_name=`hostname`
echo $host_name

if [[ "$host_name" == *borexino* ]] ; then
	CXX="g++5.1.0"; echo $CXX
else
	CXX="g++"; echo $CXX
fi

echo "Compiling $source_file"
$CXX -Wall -O3 -std=c++11 \
`root-config --cflags --glibs` \
-lEG \
-o $base_name \
$source_file