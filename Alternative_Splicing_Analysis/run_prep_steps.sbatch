#!/usr/bin/env bash
#SBATCH -A dm1_vs_sgDMPK
#SBATCH --cpus-per-task=8
#SBATCH -p fast
#SBATCH -o slurm.%j.out

module load conda
module load rmats/4.1.1
rmats_dir="/software/miniconda/envs/rmats-4.1.1/bin"  # module location

cores=8
#gtf="/shared/home/mkondili/genomes/Annotation/hg19/Homo_sapiens.GRCh37.87.CHR.gtf"
gtf="/genomes/Annotation/hg38/gencode.v39.annotation.gtf"
workDir="/dm1_vs_sgDMPK/rMATS_Splicing/"
outDir_prep="${workDir}/Prep_step/"

mkdir -p ${outDir_prep}

group1="${workDir}/Inputs/CTsgNT2_bam.txt"
group2="${workDir}/Inputs/DMsgNT2_bam.txt"
group3="${workDir}/Inputs/DMsgDMPK2_bam.txt"


###--- PREP CTsgNT2_vs_DMsgNT2 ----###
mkdir -p ${outDir_prep}/tmp_sgNT2/
python ${rmats_dir}/rmats.py --b1 ${group1} --b2  ${group2} \
--gtf ${gtf} --readLength 150 --variable-read-length \
--od ${outDir_prep} --tmp ${outDir_prep}/tmp_NTs/ \
-t "paired"  --libType fr-firststrand \
--nthread ${cores} --tstat ${cores} \
--task prep


###--- PREP DMsgNT2-vs-DMsgDMPK2 ----###
mkdir -p ${outDir_prep}/tmp_DMsgDMPK2/
python ${rmats_dir}/rmats.py --b1 ${group2} --b2  ${group3} \
--gtf ${gtf} --readLength 150 --variable-read-length \
--od ${outDir_prep} --tmp ${outDir_prep}/tmp_DMs/  \
-t "paired"  --libType fr-firststrand \
--nthread ${cores} --tstat ${cores} \
--task prep


###--- PREP CTsgNT2-vs-DMsgDMPK2 ----###
mkdir -p ${outDir_prep}/tmp_CT_vs_sgDMPK2/
python ${rmats_dir}/rmats.py --b1 ${group1} --b2  ${group3} \
--gtf ${gtf} --readLength 150 --variable-read-length \
--od ${outDir_prep} --tmp ${outDir_prep}/tmp_CT_sg2/  \
-t "paired"  --libType fr-firststrand \
--nthread ${cores} --tstat ${cores} \
--task prep
