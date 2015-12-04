#!/bin/bash
#$ -clear
#$ -S /bin/sh
#$ -w w
#$ -cwd
#$ -q 4c.q
#$ -j y
# get rid of this #$ -v PATH=$path


export LD_LIBRARY_PATH=/home/sensey/Projects/WurstImbiss/vendor/zlog/lib:$LD_LIBRARY_PATH
#valgrind --leak-check=full 
./wurstimbiss.x -l 9 -l 15 -n 50

