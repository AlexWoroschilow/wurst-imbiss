#!/bin/sh
TIMESTAMP=$(date +"%s")
ROOT="/home/stud2013/ovoroshylov/Clustering/wurst-imbiss"
export LD_LIBRARY_PATH=${ROOT}/vendor/zlog/lib:$LD_LIBRARY_PATH
${ROOT}/wurstimbiss.x -c ${ROOT}/wurstimbiss_salami.cnf > out/wurstimibiss_salami_${TIMESTAMP}.csv

