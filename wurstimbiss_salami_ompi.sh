#!/bin/bash
#0pe openmpi_pe 4
ROOT="/home/stud2013/ovoroshylov/Clustering/wurst-imbiss"
qsub -S /bin/bash -wd ${ROOT} -q "8c.q,16c.q,32c.q" ${ROOT}/wurstimbiss_salami.sh


