#!/bin/bash

source code/custom-bashrc

INPUT=$@

for element in ${INPUT[@]}
do
	echo "###### START ######" > ${element}.out
done
