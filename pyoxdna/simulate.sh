#!/bin/sh
# usage ./simulate.sh input_file

export PATH=$PATH:$OXDNA_HOME/build/bin/

trap "exit" INT TERM
trap 'kill -9 $pid > /dev/null 2>&1' EXIT # kill background process

oxDNA $1 & # run oxDNA in background
pid=$!
# echo './simulate starting '$pid
wait $! # wait for process to complete
exitCode=$?
[ $exitCode -eq 0 ] || exit 1