#!/usr/bin/env bash
#SBATCH -A dm1_vs_sgDMPK
#SBATCH --cpus-per-task=8
#SBATCH -p fast
#SBATCH -o slurm.%j.out

module load conda
module load rmats/4.1.1

rmats_dir="/software/miniconda/envs/rmats-4.1.1/bin"
cores=8
gtf="/genomes/Annotation/hg38/gencode.v39.annotation.gtf"
workDir="/dm1_vs_sgDMPK/Alternative_Splicing"

outDir_prep="${workDir}/Prep_step"
outDir_post="${workDir}/Post_step"


mkdir -p ${outDir_prep}/combo_tmp/
mkdir -p ${outDir_post}

group123="${workDir}/Inputs/all_conditions_bam.txt"


## Add a prefix to files and add them to a common tmp folder for POST task.
#> call cp_with_prefix.py
python ${rmats_dir}/../rMATS/cp_with_prefix.py "CT_vs_DMsgNT_"   ${outDir_prep}/combo_tmp/  ${outDir_prep}/tmp_NTs/*.rmats
python ${rmats_dir}/../rMATS/cp_with_prefix.py "DMsgNT_vs_DMsgDMPK2_" ${outDir_prep}/combo_tmp/  ${outDir_prep}/tmp_DMs/*.rmats
python ${rmats_dir}/../rMATS/cp_with_prefix.py "CT_vs_DMsgDMPK2_" ${outDir_prep}/combo_tmp/ ${outDir_prep}/tmp_CT_sg2/*.rmats

# Delete  in /combo_tmp/ : one set of _0, _1, _2,_3 because refer to DMsgNT samples that are already analysed first comparison
rm ${outDir_prep}/combo_tmp/CT_vs_DMsgNT_*_{0,1,2,3}.rmats #--> CTRL samples
rm ${outDir_prep}/combo_tmp/DMsgNT_vs_DMsgDMPK2_*_{0,1,2,3}.rmats #--> DM.NT
rm ${outDir_prep}/combo_tmp/CT_vs_DMsgDMPK2_*_{4,5,6,7}.rmats # --> DMsgDMPK


###----POST-treatment --------###
python ${rmats_dir}/rmats.py --b1 ${group123}  \
        --gtf ${gtf} --readLength 150 --variable-read-length  \
        --od ${outDir_post} --tmp ${outDir_prep}/combo_tmp/  \
        -t "paired"  --libType fr-firststrand   \
        --nthread ${cores} --tstat ${cores} \
        --task post


###------ STATS ------###

stats_dir=${outDir_post}/Stats_step/

##> 1/
mkdir -p ${stats_dir}/CT_vs_DMsgNT2
python ${rmats_dir}/../rMATS/rMATS_P/prepare_stat_inputs.py --new-output-dir ${stats_dir}/CT_vs_DMsgNT2 \
                    --old-output-dir ${outDir_post} --group-1-indices 0,1,2,3 --group-2-indices 4,5,6,7

python ${rmats_dir}/rmats.py --od ${stats_dir}/CT_vs_DMsgNT2 --tmp ${stats_dir}/CT_vs_DMsgNT2/tmp --task stat


##> 2/
mkdir -p ${stats_dir}/DMsgNT2_vs_DMsgDMPK2
python ${rmats_dir}/../rMATS/rMATS_P/prepare_stat_inputs.py --new-output-dir ${stats_dir}/DMsgNT2_vs_DMsgDMPK2 \
        --old-output-dir  ${outDir_post}  --group-1-indices 4,5,6,7 --group-2-indices 8,9,10,11

python ${rmats_dir}/rmats.py --od ${stats_dir}/DMsgNT2_vs_DMsgDMPK2 --tmp ${stats_dir}/DMsgNT2_vs_DMsgDMPK2/tmp  --task stat


##> 3/
mkdir -p ${stats_dir}/CT_vs_DMsgDMPK2
python ${rmats_dir}/../rMATS/rMATS_P/prepare_stat_inputs.py --new-output-dir ${stats_dir}/CT_vs_DMsgDMPK2 \
         --old-output-dir ${outDir_post} --group-1-indices 0,1,2,3 --group-2-indices 8,9,10,11

python ${rmats_dir}/rmats.py --od ${stats_dir}/CT_vs_DMsgDMPK2 --tmp ${stats_dir}/CT_vs_DMsgDMPK2/tmp --task stat
