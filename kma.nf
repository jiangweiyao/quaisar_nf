#!/usr/bin/env nextflow

Channel.fromFilePairs( params.in ).into { fastq_files; fastq_files2; fastq_files3 }
abr_ref = file(params.abrdb)
adapters = file(params.adapters)
phix = file(params.phix)
params.thread = 1

process kma_index {

    //errorStrategy 'ignore'
    publishDir params.out, overwrite: true

    output:
    path "abr*" into kma_index_out

    """
    kma index -i ${abr_ref} -o abr
    """
}

//kma_index_out.print()

process kma_map {

    //errorStrategy 'ignore'
    publishDir params.out, mode: 'copy', overwrite: true

    input:
    path index from kma_index_out
    tuple val(name), file(fastq) from fastq_files2

    output:
    path "*_abr*" into kma_out

    //kma -i ${fastq[0]} -ipe ${fastq[1]} -o ${name}_abr.out -t_db abr -1t1


    """
    kma -ipe ${fastq} -o ${name}_abr -t_db abr -1t1
    """
}

