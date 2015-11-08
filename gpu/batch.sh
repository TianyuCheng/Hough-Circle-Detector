#!/bin/tcsh

SBATCH -J hcd              # Job Name
SBATCH -o hcd.o%j          # Output and error file name
SBATCH -n 1                # Total number of GPUs requested
SBATCH -p gpudev           # Queue name
SBATCH -t 00:10:00         # Run time (hh:mm:ss)
SBATCH -A tcheng

./main
