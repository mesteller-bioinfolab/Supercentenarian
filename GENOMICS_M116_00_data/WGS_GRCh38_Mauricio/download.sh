#!/bin/bash
#
#SBATCH -p normal # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J download
#SBATCH --mem 80G # memory pool for all cores
#SBATCH -t 1-00:00 # time (D-HH:MM)
#SBATCH -o ./log/log.%j.out # STDOUT
#SBATCH -e ./log/log.%j.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=eva.brigos@upf.edu # send-to address

wget ftp://FTPuser:FTPusersPassword@xfer13.crg.eu:221/external_files/supercentenarian/*

