#!/bin/bash

[[ $# == 1 ]] && ( [ "$1" == "egst" ] || [ "$1" == "gst" ] || [ "$1" == "gtrac" ] ) || ( echo "$0 [egst|gst|gtrac]" && exit )

base_name="$1Tally"
source_file="$1Tally.cc"
host_name=`hostname`
echo $host_name

if [[ "$host_name" == *borexino* ]] ; then
	CXX="g++6.2.0 -Wl,-rpath,/storage/gpfs_data/borexino/users/viatro/opt/gcc-6.2.0/lib64"; echo $CXX
else
	CXX="g++"; echo $CXX
fi

echo "Compiling $source_file"
$CXX -Wall -O3 -std=c++14 \
`root-config --cflags --glibs` \
-lEG \
-o $base_name \
$source_file