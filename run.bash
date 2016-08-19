#!/bin/bash

source code/custom-bashrc


if [ "$#" -ne 2 ]; then
        printf "\n\n###### Usage\n\n"
        printf "$0 <PATH to params file> <PATH to group_info.txt file>\n\n"
        exit
fi

PARAMS=$1
GROUP_FILE=$2

samples=`cut -f1 $GROUP_FILE | code/skipn 1`
code/setup_waiting.bash ${samples[*]}
echo ${samples[*]} | xargs -n1 code/job_submitter.bash 5 code/By_Sample $PARAMS
code/wait_jobs.bash ${samples[*]}


GROUP_BYs=`head -1 $GROUP_FILE | cut -f2-`
code/setup_waiting.bash ${GROUP_BYs[*]}
echo ${GROUP_BYs[*]} | xargs -n1 code/job_submitter.bash 6 code/Summarize $PARAMS pipeline/summarize $GROUP_FILE
code/wait_jobs.bash ${GROUP_BYs[*]}

