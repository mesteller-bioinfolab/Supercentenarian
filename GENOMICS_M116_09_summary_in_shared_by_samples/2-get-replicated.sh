#!/bin/bash

rm -rf ./real/tmp/*

COMMAND="real/scripts/get-replicated.sh"


sbatch ${COMMAND}

