#!/bin/bash
# Copyright 2015 Alex Woroschilow (alex.woroschilow@gmail.com)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#$ -clear
#$ -q stud.q
#$ -S /bin/bash
ROOT=$(pwd)
TIMESTAMP=$(date +"%s")
echo "working folder: ${ROOT}"
export LD_LIBRARY_PATH=${ROOT}/vendor/zlog/lib:$LD_LIBRARY_PATH
${ROOT}/wurstimbiss.x -c ${ROOT}/etc/structure.cnf > ${ROOT}/out/structure_${TIMESTAMP}.csv