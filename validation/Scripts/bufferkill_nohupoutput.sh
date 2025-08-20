#!/bin/bash

## This script is meant to look for active processes with g4cmpquasiparticle as the process name,
## and if it finds any, periodically kill the text in the nohup.out file without killing the file
## This is meant to allow us to look at the output of a run without having to actively monitor
## and kill that text ourselves.
## This needs to be executed in the same directory as that in which the nohup file is located

#Get the number of processes running with g4cmpquasiparticle in the name.
QPPROCESSES=`ps | grep g4cmpQuasiparticle | wc -l`
while true
do
    sleep 1
    if [ $QPPROCESSES -gt 1 ]
    then
	sleep 5
	SIZE=`stat -f "%z" nohup.out`
	sleep 0.5
	SIZE2=`stat -f "%z" nohup.out`
	if [ $SIZE -ne $SIZE2 ]
	then
	    > nohup.out
	fi
    fi
done
  
