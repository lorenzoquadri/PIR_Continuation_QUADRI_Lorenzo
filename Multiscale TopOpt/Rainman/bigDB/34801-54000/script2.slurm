#!/bin/bash

#SBATCH --nodes=4
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --time=48:00:00
#SBATCH --begin=now
#SBATCH --mail-user=lorenzo.quadri@student.isae-supaero.fr
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=PIR2_LQ
#SBATCH -o slurm.%j.out # STDOUT
#SBATCH -e slurm.%j.err # STDERR

module load matlab/2018a

matlab -r "jojo=1;run('generateDataFromGrid.m');exit;"