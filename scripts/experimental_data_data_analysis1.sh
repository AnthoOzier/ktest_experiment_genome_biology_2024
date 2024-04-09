#!/bin/bash  
#SBATCH --exclusive   
#SBATCH --mail-user=____@___.__ 
#SBATCH --mail-type=all,abort,end,time_limit 
#SBATCH -o job_ccdfK-%j.out 
#SBATCH --export=ALL
#  #SBATCH --nodes=2



module purge
module load cuda/9.2.148
module load gcc/10.1
module load cmake/3.10.2
source miniconda3/bin/activate
conda activate mmd3

path=$1
diroutput=$2
dirlog=$3
jn=$4
int=$5
bool=$6
str=$7 
nj=$8
partition=$9

echo ${path}

# if [ ! -e ${diroutput}${jn}.csv ]; then 
echo "${jn} launched with"$'\n'"str=${str}"$'\n'"bool=${bool}"$'\n'"int=${int} parallel : ${nj}"
srun -n 1 ~/miniconda3/envs/mmd3/bin/python ${path}experimental_data_data_analysis.py -s $str -b $bool -i $int >${dirlog}log1_${jn} 2>&1
# fi  
