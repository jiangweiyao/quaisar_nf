#!/usr/bin/env nextflow

Channel.fromFilePairs( params.in ).into { fastq_files; fastq_files2; fastq_files3; fastq_files4; fastq_files5 }
abr_ref = file(params.abrdb)
adapters = file(params.adapters)
phix = file(params.phix)
params.thread = 1


mash_genome_db = file(params.genome_db)
mash_plasmid_db = file(params.plasmid_db)

process mash_screen_genome {

    errorStrategy 'ignore'
    publishDir params.out, overwrite: true
    maxForks 1

    input:
    tuple val(name), file(fastq) from fastq_files4

    output:
    path "*_pathogen_id.out" into mash_screen_genome_out

    """
    cat ${fastq[0]} ${fastq[1]} > combined.fastq.gz
    mash screen -w ${mash_genome_db} combined.fastq.gz | sort -gr > ${name}_pathogen_id.out
    """
}

//mash_screen_genome_out.print()

/*
process mash_screen_plasmid {

    errorStrategy 'ignore'
    publishDir params.out, overwrite: true

    input:
    tuple val(name), file(fastq) from fastq_files5

    output:
    path "*_plasmid_id.out" into mash_screen_plasmid_out

    """
    cat ${fastq[0]} ${fastq[1]} > combined.fastq.gz
    mash screen -w ${mash_plasmid_db} combined.fastq.gz | sort -gr > ${name}_pathogen_id.out
    """
}*/


