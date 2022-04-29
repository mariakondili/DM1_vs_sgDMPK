#!/usr/bin/env nextflow

// ~~~ Subject ~~
// Align sequences of human cells in Control, DM1-w/o treatment & DM1-with-treatment-sg2DMPK conditions,
// Count gene-expression and splicing events

// Author: Maria Kondili
// Date : 24 January 2022

workflow.onComplete {
    //any worlflow property can be used here
    if ( workflow.success ) {
        println "Pipeline Complete"
    }
    println "Command line: $workflow.commandLine"
}

workflow.onError {
    println "Oops .. something went wrong"
}


params.samplelist   = "/dm1_vs_sgDMPK/Nextflow_Pip/test_samplesheet.tsv"
params.outDir       = "/dm1_vs_sgDMPK/BAM"
params.data_dir     = "/dm1_vs_sgDMPK/Fastq/raw/Merged_Lanes/"
params.cutadapt_dir = "/dm1_vs_sgDMPK/Fastq/CutAdapted/"
params.countsDir    = "/dm1_vs_sgDMPK/HTSeq_counts/"
params.refGenome    = "/shared/bank/homo_sapiens/hg38/star-2.7.5a"
params.ref_Gtf      = "/shared/genomes/Annotation/hg38/gencode.v39.annotation.gtf"
params.cpus = 20

// create pairs as variables
samplist = file("${params.samplelist}")
reader = samplist.newReader()
// read the samplesheet and create array with two reads to Input in next Process
pair=[]
samplist.withReader {
    String line
    while ( line = reader.readLine() ) {
        String r1 = line.split(" ")[0]  // press Tab to obtain space.
        println "read_1 = ${r1}"
        String r2 = line.split(" ")[1]
        println "read_2 = ${r2}"
        String bn = r1.split("${params.data_dir}")[1].split("_R1.fastq.gz")[0] //.split("_R1_001")[0]
        println "Base-name = ${bn}"
        pair.add([bn,r1,r2])
    }
}

inputChannel  = Channel.fromList(pair)

process Cut_Adapters {

        cpus "${params.cpus}"
        memory "15G"
        module "cutadapt/2.10"

        input:
        tuple val(bn),path(r1),path(r2) from inputChannel

        output:
        tuple val(bn),path("${bn}_R1_cutadapt.fastq"),path("${bn}_R2_cutadapt.fastq") into cutadChannel

        shell:
        """
        module li ;
        mkdir -p !{params.cutadapt_dir}/!{bn}/

        cutadapt --cores !{params.cpus} \
        -B CTGTCTCTTATA \
        --minimum-length 75 --error-rate 0.1 \
        --report  minimal \
        --output !{bn}_R1_cutadapt.fastq \
        -p       !{bn}_R2_cutadapt.fastq \
	      !{r1} !{r2}

        cp !{bn}_R1_cutadapt.fastq !{params.cutadapt_dir}/!{bn}/
        cp !{bn}_R2_cutadapt.fastq !{params.cutadapt_dir}/!{bn}/

        """
}



process STAR_Alignment {

      cpus "${params.cpus}"
      memory "40G"
      module "star/2.7.5a:perl/5.26.2"

      input:
      tuple val(bn),path(f1),path(f2) from cutadChannel

      output:
      tuple path("${bn}_Aligned.sortedByCoord.out.bam"),val(bn) into alignChannel

      shell:
      """
      mkdir -p !{params.outDir};
      echo "SampleID=" !{bn};
      STAR --runMode  alignReads \
           --runThreadN !{params.cpus} \
           --genomeDir !{params.refGenome} \
           --readFilesIn !{f1} !{f2}  \
           --outFileNamePrefix !{bn}_  \
           --outSAMtype BAM SortedByCoordinate \
           --outReadsUnmapped Fastx ;
           ## --readFilesCommand zcat | gunzip -c -> use unzipped fastq
	   cp !{bn}_Aligned.sortedByCoord.out.bam  !{params.outDir}/ ;
	   cp !{bn}_Log.final.out !{params.outDir}/

       """
}


process HTSeq_counts {

      cpus "${params.cpus}"
    	memory "20G"
    	module "samtools/1.10:htseq/0.12.4"

    	input:
    	tuple path(bam),val(bn) from alignChannel

    	output:
    	path("${bn}_counts.txt") into htseqChannel

    	shell:
        """
        mkdir -p !{params.countsDir};

        samtools index !{bam} ;

        if [ -e !{bam}.bai ] ; then echo "Index created"; fi

        htseq-count -n !{params.cpus} \
        -r pos -s reverse \
        --type gene \
        --idattr gene_id \
        --mode union --minaqual 10 \
        --counts_output !{bn}_counts.txt \
        !{bam}  !{params.ref_Gtf} ;

        cp  !{bn}_counts.txt  !{params.countsDir}/

        """
}


// LAUNCH PIPE:
// nextflow -c  nextflow.confing  run -w  NXF_pip/  main.nf
