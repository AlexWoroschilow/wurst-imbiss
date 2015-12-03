#!/bin/bash
#$ -clear
#$ -S /bin/sh
#$ -w w
#$ -cwd
#$ -q 4c.q
#$ -j y
# get rid of this #$ -v PATH=$path



#valgrind --leak-check=full 
./wurstimbiss.x -a abc/alphabet_t500 -s seq/pdb_seq.list -b ./pdb/pdb_02.lib -l 9 -l 15 -n 50

