#!/usr/bin/env nextflow


Channel.fromFilePairs( params.in ).into { fastq_files; fastq_files2; fastq_files3 }
abr_ref = file(params.abrdb)
adapters = file(params.adapters)
phix = file(params.phix)
params.thread = 1


process fastqc {
    
    errorStrategy 'ignore'
    publishDir params.out, pattern: "*.html", mode: 'copy', overwrite: true

    input:
    set val(name), file(fastq) from fastq_files
 
    output:
    file "*_fastqc.{zip,html}" into qc_files, qc_files1

    """
    fastqc -q ${fastq}
    """
}

//qc_files1.collect().print()

process multiqc {

    errorStrategy 'ignore'
    publishDir params.out, mode: 'copy', overwrite: true

    input:
    file reports from qc_files.collect().ifEmpty([])

    output:
    path "multiqc_report.html" into multiqc_output

    """
    multiqc $reports
    """
}

process srst2 {

    errorStrategy 'ignore'
    //maxForks 1
    publishDir params.out, mode: 'copy', overwrite: true

    input:
    tuple val(name), file(fastq) from fastq_files2

    output:
    path "*results.txt" into srst2_output

    """
    srst2 --input_pe ${fastq} --output ${name}_srst2 --log --gene_db ${abr_ref}
    """
}

process adapter_trimming {

    //errorStrategy 'ignore'
    //publishDir params.out, mode: 'copy', overwrite: true

    input:
    tuple val(name), file(fastq) from fastq_files3

    output:
    tuple val(name), file("*_clean{1,2}.fq.gz") into trimmed_fastq

    """
    bbduk.sh -Xmx1g in1=${fastq[0]} in2=${fastq[1]} out1=int1.fq.gz out2=int2.fq.gz ref=${phix} k=31 hdist=1 
    bbduk.sh -Xmx1g in1=int1.fq.gz in2=int2.fq.gz out1=${name}_clean1.fq.gz out2=${name}_clean2.fq.gz ref=${adapters} ktrim=r k=23 mink=11 hdist=1 tpe tbo
    """
}

process assembly {

    //errorStrategy 'ignore'
    //publishDir params.out, mode: 'copy', overwrite: true

    input:
    tuple val(name), file(fastq) from trimmed_fastq

    output:
    tuple val(name), path("*_scaffolds.fasta") into assembly_output

    """
    spades.py -1 ${fastq[0]} -2 ${fastq[1]} -o ${name} -m 6 -t ${params.thread}
    cp ${name}/scaffolds.fasta ${name}_scaffolds.fasta
    """
}

process assembly_size_filter {

    errorStrategy 'ignore'
    publishDir params.out, mode: 'copy', overwrite: true

    input:
    tuple val(name), file(assembly) from assembly_output

    output:
    tuple val(name), path("*_scaffolds_filtered.fasta") into assembly_filter_output, assembly_filter_output2, assembly_filter_output3

    """
    reformat.sh in=${assembly} out=${assembly.simpleName}_filtered.fasta minlength=${params.sizefilter}
    """
}

process prokka {

    errorStrategy 'ignore'
    publishDir params.out, mode: 'copy', overwrite: true

    input:
    tuple val(name), file(assembly) from assembly_filter_output

    output:
    tuple val(name), path("*") into prokka_output

    """
    prokka --cpus ${params.thread} --outdir ${name}_prokka --prefix ${name} ${assembly}
    """
}

//assembly_filter_output2.print()

process quast {

    errorStrategy 'ignore'
    publishDir params.out, mode: 'copy', overwrite: true

    input:
    tuple val(name), file(assembly) from assembly_filter_output2

    output:
    path("*") into quast_output

    """
    quast ${assembly} -o ${name}_quast
    """
}

process busco {

    //errorStrategy 'ignore'
    publishDir params.out, mode: 'copy', overwrite: true

    input:
    tuple val(name), file(assembly) from assembly_filter_output3

    output:
    path("*") into busco_output

    """
    busco --auto-lineage-prok -f -m geno -o ${name}_busco -i ${assembly}
    """
}

