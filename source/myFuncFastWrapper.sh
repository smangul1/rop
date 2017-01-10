#!/bin/bash
# myFuncFastWrapper.sh
echo $SGE_TASK_ID
./myFunc.sh $SGE_TASK_ID
