#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

#$ -N zacros_JA 					#This is the name of the job array
#$ -t 1-10  							#Assumes task IDs increment by 1; can also increment by another value
#$ -tc 10 							#This is the total number of tasks to run at any given moment
#$ -pe openmpi-smp 2 -l mem_free=1G			#Change the last field to the number of processors desired per task

job_file='/mnt/files-ccei/home/wangyf/Ceria/Pd_Growth/3d/all/61/300/1/ss/4/outputs/dir_list.txt'
#Change to the job directory
job_path=$(sed -n "$SGE_TASK_ID p" "$job_file")
cd "$job_path" #SGE_TASK_ID is the task number in the range <task_start_index> to <task_stop_index>


time /home/vlachos/wangyf/programs/zacros_2.1/build/zacros.x