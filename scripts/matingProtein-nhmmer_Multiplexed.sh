#!/bin/bash

#SBATCH -J matingProteinNhmmer #Job name for the array
#SBATCH -o Diagnostic/matingProteinNhmmerMultiplexed_%A-%a.out #File to which standard out will be written
#SBATCH -e Diagnostic/matingProteinNhmmerMultiplexed_%A-%a.err #File to which standard err will be written
#SBATCH -p serial_requeue #Partition to submit to
#SBATCH -t 22:30:00 #Runtime in HH:MM:SS
#SBATCH --mem-per-cpu=10240 #Memory per cpu in MB (see also --mem)
#SBATCH -n 1 #Number of cores

#export SLURM_ARRAY_TASK_ID=test
#echo "$SLURM_ARRAY_TASK_ID.fasta"
##To multiplex accross chromosomes of all the genomes, I am going to build a fileList and submit arrayed jobs with the appropriate file name from the list.

varCount=0
echo Searching genomes from $1
#for file in "$1"/*
while read -r line
do
        ##I need to make an array of filenames.
        #((index=varCount - 1))
        #echo $line
        fileList[$varCount]="$line"
        ((varCount=varCount + 1))
done < $2

#echo $((${SLURM_ARRAY_TASK_ID} + 1))
#echo $((10000*$4))
index=$((10000*${4} + ${SLURM_ARRAY_TASK_ID}))
echo Working on genome saved in file numbered $index
echo Working on ${fileList[${index}]} from batch ${4}
nhmmer --qfasta --tblout nhmmerOutput/tbl-nhmmer_${3}_${fileList[${index}]} -o nhmmerOutput/nhmmer_${3}_${fileList[${index}]} $3 $1/${fileList[${index}]}

#$SLURM_ARRAY_TASK_ID sets the folder containing the genomes as FASTA files to be searched in by nhmmer using the fileList stored before script run in $2. $3 is the FASTA file containing the DNA sequence(s) that is to be found in the genomes; and the outputs' names will be a combination of the genes to be searched for and the genomes searched in. $4 is the job multiplier since SLURM has a MaxArraySize on Odyssey of 10000, the true file index is therefore ((10000*$2 + ${SLURM_ARRAY_TASK_ID})).
#Use "$sbatch --array=<start-end> <script-name> to run. Good luck!

###SBATCH --mail-type=END #Type of email notification- BEGIN,END,FAIL,ALL
###SBATCH --mail-user=sriramsrikant@fas.harvard.edu #Email to which notifications will be sent
