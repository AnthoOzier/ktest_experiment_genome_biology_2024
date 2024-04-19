
path=/home/scripts/ 
dirdata=/home/data/experimental_data/
n_jobs=64
t=01:00:00
pzstr=''
data_str=2023_06_ccdf${pzstr}
stat="kfda"
nystrom=false
kernel="gauss" #"fisher_zero_inflated_gaussian" "naive_fisher_zero_inflated_gaussian" "linear" "gauss_quantile"
quantile=80 # if kernel == "gauss_quantile"
ns="20 40 60 80 100 160 200"
permutation=false
ignore_zeros=false
n_permutations=1000
lots=1000

if [ "$kernel" = "gauss" ]; then
    kernel_str=''
elif [ "$kernel" = "linear" ]; then
    kernel_str='_klin_'
elif [ "$kernel" = "fisher_zero_inflated_gaussian" ]; then
    kernel_str='_kfzig_'
elif [ "$kernel" = "gauss_quantile" ]; then
    kernel_str=_kgq${quantile}_
else
    kernel_str='_naivekfzig_'
fi

if [ "$permutation" = "true" ]; then
    perm_str='_perm_'
else
    perm_str=''
fi

if [ "$nystrom" = "true" ]; then
    ny_str='_ny_'
else
    ny_str=''
fi
if [ "$ignore_zeros" = "true" ]; then
    ignore_zeros_str='_iz_'
else
    ignore_zeros_str=''
fi

model_str=${stat}${kernel_str}${ny_str}${perm_str}${ignore_zeros_str}

diroutput=/home/data/experimental_data/${data_str}_${model_str}/
dirlog=/home/data/experimental_data/log_${data_str}_${model_str}/
dirscores=/home/data/experimental_data/scores_${data_str}_${model_str}/

mkdir -p ${diroutput} 
mkdir -p ${dirlog} 
mkdir -p ${dirscores} 

for n in $ns; do       
    echo n${n} t${t//:/_} ${data_str} ${model_str}
    nseeds=10
    seeds=$(seq 0 $nseeds 90)
    hundreds="0 100 200 300 400"
    for seed in $seeds; do
        for hundred in $hundreds; do  
            bool=nystrom:${nystrom}+permutation:${permutation}+ignore_zeros:${ignore_zeros}
            s0=$((hundred+seed))
            s1=$((hundred+seed+nseeds))
            int=n:${n}+n_jobs:${n_jobs}+s0:${s0}+s1:${s1}+lots:${lots}+n_permutations:${n_permutations}+quantile:${quantile}
            jobname=n${n}s${s0}-${s1}_${data_str}_${model_str}_lots${lots}_t${t//:/_}_njobs${n_jobs}_${partition}
            filename=${stat}_n${n}
            str=stat:${stat}+diroutput:${diroutput}+dirdata:${dirdata}+filename:${filename}+kernel:${kernel}
            echo run ${jobname}
            sbatch --job-name=${jobname} --ntasks-per-node=$n_jobs --time=$t -p $partition ${path}experimental_data_data_analysis1.sh $path $diroutput $dirlog $jobname $int $bool $str $njobs $partition >${dirlog}log0_${jobname} 2>&1
        done
    done
done


