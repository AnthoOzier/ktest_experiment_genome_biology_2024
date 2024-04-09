

path=/home/scripts/

files="Angelidis2019_pneumo Angelidis2019_alvmac Reyfman2020_alvmac Reyfman2020_pneumo"
#files="CanoGamez2020_Memory-iTreg CanoGamez2020_Memory-Th0 CanoGamez2020_Memory-Th17 CanoGamez2020_Memory-Th2 CanoGamez2020_Naive-iTreg CanoGamez2020_Naive-Th0 CanoGamez2020_Naive-Th17 CanoGamez2020_Naive-Th2 Hagai2018_pig-lps Hagai2018_mouse-lps Hagai2018_mouse-pic Hagai2018_rat-pic Hagai2018_rabbit-lps Hagai2018_rat-lps"
n_jobs=64
ignore_zeros=false
binary=false
kernel="gauss" #"fisher_zero_inflated_gaussian" "naive_fisher_zero_inflated_gaussian" "linear" "gauss_quantile"
quantile=80 # if kernel == "gauss_quantile"
sn=true # analysis on the normalized data 
center_non_zero_only=true # correct for the batch effect on the non-zero observations only
center_by="#-replicate_#plus#label_in_input_space" # Adapted to Angelidis and Reyfman datasets
center_by='replicate_in_input_space' # Adapted to others
multivariate=false
lots=500   
nvariables=1000

if [ "$kernel" = "gauss" ]; then
    kernel_str=''
elif [ "$kernel" = "linear" ]; then
    kernel_str='klin_'
elif [ "$kernel" = "gauss_quantile" ]; then
    kernel_str=_kgq${quantile}_
else
    kernel_str='kfzig_'
fi

if [ "$sn" = "true" ]; then
    sn_str=''
else
    sn_str='sct_'
fi

if [ "$multivariate" = "true" ]; then
    mv_str='mv'
else
    mv_str=''
fi

if [ "$ignore_zeros" = "true" ]; then
    iz_str='iz_'
else
    iz_str=''
fi

if [ "$center_non_zero_only" = "true" ]; then
    cbnz_str='cbnz_'
else
    cbnz_str=''
fi

if [ "$center_by" = "replicate" ]; then
    cb_str='cbr_'
elif [ "$center_by" = "#-replicate_#plus#label" ]; then
    cb_str='cbrl_'
elif [ "$center_by" = "replicate_in_input_space" ]; then
    cb_str='cbris_'
elif [ "$center_by" = "#-replicate_#plus#label_in_input_space" ]; then
    cb_str='cbrlis_'
else
    cb_str=''
fi

model_str=${sn_str}${mv_str}${iz_str}${cb_str}${cbnz_str}${kernel_str}
diroutput=/home/data/squair/Courtine_KFDA_${model_str}results/
dirlog=/home/data/squair/Courtine_KFDA_${model_str}log/
dirdata=/home/data/squair/normalized/

mkdir -p ${diroutput} # créer le dossier des outputs
mkdir -p ${dirlog} # créer le dossier des outputs



for file in $files; do
    if [ "$sn" = "true" ]; then 
        file=normalized_${file}
    else
        file=residuals_${file}
    fi

    variables=$(seq 0 $nvariables 19100) 

    if [ "$file" = "normalized_Angelidis2019_pneumo" ]; then
        t=00:15:00
        variables=$(seq 0 $nvariables 16300)
    fi
    if [ "$file" = "normalized_Angelidis2019_alvmac" ]; then
        t=00:20:00
        variables=$(seq 0 $nvariables 14500)
    fi

    if [ "$file" = "normalized_Reyfman2020_pneumo" ]; then
        n_jobs=3
        t=06:40:00
        variables=$(seq 0 $nvariables 25600)
    fi
    if [ "$file" = "normalized_Reyfman2020_alvmac" ]; then
        n_jobs=3
        t=04:40:00
        variables=$(seq 0 $nvariables 24600)
    fi
    if [ "$file" = "normalized_Hagai2018_pig-lps" ]; then
        t=01:00:00
        variables=$(seq 0 $nvariables 8900)
    fi
    if [ "$file" = "normalized_Hagai2018_rabbit-lps" ]; then
        t=03:00:00
        n_jobs=10
        variables=$(seq 0 $nvariables 9200)
    fi
    if [ "$file" = "normalized_Hagai2018_rat-lps" ]; then
        t=00:50:00
        variables=$(seq 0 $nvariables 14900)
    fi
    if [ "$file" = "normalized_Hagai2018_rat-pic" ]; then
        t=00:50:00
        variables=$(seq 0 $nvariables 14900)
    fi
    if [ "$file" = "normalized_Hagai2018_mouse-pic" ]; then
        t=02:00:00
        variables=$(seq 0 $nvariables 16000)
    fi
    if [ "$file" = "normalized_Hagai2018_mouse-lps" ]; then
        n_jobs=20
        variables=$(seq 0 $nvariables 15000)
        t=03:00:00  
    fi   
    if [ "$file" = "normalized_CanoGamez2020_Memory-Th0" ]; then
        variables=$(seq 0 $nvariables 19100)
        t=00:40:00
    fi   
    if [ "$file" = "normalized_CanoGamez2020_Memory-Th2" ]; then
        t=00:40:00
        variables=$(seq 0 $nvariables 18100)
    fi   
    if [ "$file" = "normalized_CanoGamez2020_Memory-Th17" ]; then
        t=00:40:00
        variables=$(seq 0 $nvariables 18100)
    fi   
    if [ "$file" = "normalized_CanoGamez2020_Naive-Th0" ]; then
        t=00:20:00
        variables=$(seq 0 $nvariables 17100)
    fi  
    if [ "$file" = "normalized_CanoGamez2020_Naive-iTreg" ]; then
        t=00:30:00
        variables=$(seq 0 $nvariables 18100)
    fi  
    if [ "$file" = "normalized_CanoGamez2020_Memory-iTreg" ]; then
        t=00:40:00
        variables=$(seq 0 $nvariables 19100)
    fi  
    if [ "$file" = "normalized_CanoGamez2020_Naive-Th2" ]; then
        t=00:30:00
        variables=$(seq 0 $nvariables 19100)
    fi   
    if [ "$file" = "normalized_CanoGamez2020_Naive-Th17" ]; then
        t=00:40:00
        variables=$(seq 0 $nvariables 19100)
    fi   

    if [ "$multivariate" = "false" ]; then
        for v in $variables; do
            v0=$v
            v1=$((v+nvariables))
            diroutputv=${diroutput}${v0}_${v1}       
            bool=ignore_zeros:${ignore_zeros}+binary:${binary}+multivariate:${multivariate}+center_non_zero_only:${center_non_zero_only}
            int=n_jobs:${n_jobs}+time:${t//:/_}+v0:${v0}+v1:${v1}+lots:${lots}+quantile:${quantile}
            jobname=${file}_${sn_str}${kernel_str}${cb_str}${cbnz_str}${iz_str}_v${v0}_${v1}_t${t//:/_}_${partition}
            filename=DEA_${file}
            str=file:${file}+diroutput:${diroutputv}+dirdata:${dirdata}+filename:${filename}+partition:${partition}+kernel:${kernel}+center_by:${center_by}
            echo run ${jobname}
            sbatch --job-name=${jobname} --time=$t -p $partition ${path}Courtine_DEA1.sh $path $diroutputv $dirlog $jobname $int $bool $str $njobs $partition 
        done
    else
        v0=0
        v1=0
        t=00:06:00
        diroutputv=${diroutput} 
        bool=ignore_zeros:${ignore_zeros}+binary:${binary}+multivariate:${multivariate}+center_non_zero_only:${center_non_zero_only}
        int=n_jobs:${n_jobs}+time:${t//:/_}+v0:${v0}+v1:${v1}+lots:${lots}
        jobname=${file}_${kernel}_mv_v${v0}_${v1}_cb${center_by}_t${t//:/_}_${partition}
        filename=DEA_${file}
        str=file:${file}+diroutput:${diroutputv}+dirdata:${dirdata}+filename:${filename}+partition:${partition}+kernel:${kernel}+center_by:${center_by}
        echo run ${jobname}
        sbatch --job-name=${jobname} --time=$t -p $partition ${path}squair_data_analysis1.sh $path $diroutputv $dirlog $jobname $int $bool $str $njobs $partition 
    fi
done



