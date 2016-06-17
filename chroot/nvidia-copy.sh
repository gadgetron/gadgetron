#! /bin/sh
# Copyright (c) 2015, NVIDIA CORPORATION. All rights reserved.
# Modified from nvidia-docker by Michael S. Hansen (michael.hansen@nih.gov

NV_DEVICE="/dev/nvidia"
UVM_DEVICE="${NV_DEVICE}-uvm"
CTL_DEVICE="${NV_DEVICE}ctl"

NV_BINS_VOLUME="/usr/local/bin"
NV_BINS="nvidia-cuda-mps-control \
         nvidia-cuda-mps-server \
         nvidia-debugdump \
         nvidia-persistenced \
         nvidia-smi"

NV_LIBS_VOLUME="/usr/local/nvidia"
NV_LIBS_CUDA="cuda \
              nvcuvid \
              nvidia-compiler \
              nvidia-encode \
              nvidia-ml"

NV_DOCKER_ARGS=""

__log()
{
    local level="$1"
    local msg="$2"

    printf "[ NVIDIA ] =$level= $msg\n" >&2
}

__library_paths()
{
    local lib="$1"

    echo $( ldconfig -p | grep "lib${lib}.so" | awk '{print $4}' )
}

__library_arch()
{
    local lib="$1"

    echo $( file -L $lib | awk '{print $3}' | cut -d- -f1 )
}

__filter_duplicate_paths()
{
    local paths="$1"

    local sums="$( md5sum $paths | sed 's/[^/]*$/ &/' )"
    local uniq="$( echo "$sums" | uniq -u -f2 | awk '{print $2$3}')"
    local dupl="$( echo "$sums" | uniq --all-repeated=separate -f2 \
                                | uniq -w 32 | awk 'NF {print $2$3}')"
    echo $uniq $dupl
}

copy_files()
{
    local dest="$1" 
    mkdir -p ${dest}/$NV_LIBS_VOLUME/lib
    mkdir -p ${dest}/$NV_LIBS_VOLUME/lib64
    for lib in $NV_LIBS_CUDA; do
        local paths="$( __library_paths $lib )"
        if [ -z "$paths" ]; then
            __log WARN "Could not find library: $lib"
            continue
        fi
        for path in $( __filter_duplicate_paths "$paths" ); do
	   echo "cp $path ${dest}$path"
	   cp $path ${dest}$path
            case $( __library_arch "$path" ) in
                32) echo "cp $path ${dest}$NV_LIBS_VOLUME/lib/$(basename $path)";cp $path ${dest}$NV_LIBS_VOLUME/lib/$(basename $path) ;;
                64) echo "cp $path ${dest}$NV_LIBS_VOLUME/lib64/$(basename $path)";cp $path ${dest}$NV_LIBS_VOLUME/lib64/$(basename $path) ;;
            esac
        done
    done

    mkdir -p ${dest}/$NV_BINS_VOLUME/lib
    for bin in $NV_BINS; do
        local path="$( which $bin )"
        if [ -z $path ]; then
            __log WARN "Could not find binary: $bin"
            continue
        fi
	echo "cp $path ${dest}${NV_BINS_VOLUME}/$bin"
	cp $path ${dest}${NV_BINS_VOLUME}/$bin
    done
}

echo "Copying files...."
copy_files $1
echo "Copying done."