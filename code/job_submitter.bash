#!/bin/bash

source code/custom-bashrc

NSLOTS=$1
shift

qsub -b Y -cwd -pe threaded $NSLOTS -e /dev/null -o /dev/null -l tmp_free=120G -l tmp_token=30G -hard $@
