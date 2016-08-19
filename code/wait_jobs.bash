#!/bin/bash

source code/custom-bashrc

INPUT=$@

TOTAL=`echo ${INPUT[*]} | sed 's/\s\+/\n/g' | wc -l`

STATUS=`echo ${INPUT[*]} | sed 's/\s\+\|$/.out /g' | xargs -n1 tail -1 | grep "^###### Done ######$" | wc -l`

while [ $STATUS -lt $TOTAL ]
do
	STATUS=`echo ${INPUT[*]} | sed 's/\s\+\|$/.out /g' | xargs -n1 tail -1 | grep "^###### Done ######$" | wc -l`
	sleep 10
done
