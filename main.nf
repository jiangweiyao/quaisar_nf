#!/usr/bin/env nextflow


Channel.fromFilePairs( params.in ).into { fastq_files; fastq_files2; fastq_files3; fastq_files4; fastq_files5  }
abr_ref = file(params.abrdb)
adapters = file(params.adapters)
phix = file(params.phix)
params.thread = 1
mash_genome_db = file(params.genome_db)
mash_plasmid_db = file(params.plasmid_db)

println """\
         Q U A I S A R - H     NEXTFLOW      PIPELINE   
         =====================================================
         input reads (--in)                  : ${params.in}
         outdir (--out)                      : ${params.out}
         Antibiotic Resistance DB (--abrdb)  : ${params.abrdb}
         adapters (--adapters)               : ${params.adapters}
         Mash Genome Reference (--genome_db) : ${params.genome_db}
         size cutoff (--sizefilter)          : ${params.sizefilter}
         """
         .stripIndent()


// Check if Mash reference file already exists. If not, download it.
mash_genome_file = file("$baseDir/db/refseq.genomes.k21s1000.msh")
if(!mash_genome_file.exists()){
    println("Mash genome reference missing. Downloading...")
    mash_genome_file = file('https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh')
    mash_genome_file.copyTo("$baseDir/db/refseq.genomes.k21s1000.msh")
}

// Check if Mash plasmid file already exists. If not, download it.
mash_plasmid_file = file("$baseDir/db/refseq.plasmid.k21s1000.msh")
if(!mash_plasmid_file.exists()){
    println("Mash plasmid reference missing. Downloading...")
    mash_plasmid_file = file('https://gembox.cbcb.umd.edu/mash/refseq.plasmid.k21s1000.msh')
    mash_plasmid_file.copyTo("$baseDir/db/refseq.plasmid.k21s1000.msh")
}

process fastqc {
    
    //errorStrategy 'ignore'
    publishDir params.out, pattern: "*.html", mode: 'copy', overwrite: true

    input:
    set val(name), file(fastq) from fastq_files
 
    output:
    file "*_fastqc.{zip,html}" into qc_files, qc_files1

    """
    fastqc -q ${fastq}
    """
}

process multiqc {

    //errorStrategy 'ignore'
    publishDir params.out, mode: 'copy', overwrite: true

    input:
    file reports from qc_files.collect().ifEmpty([])

    output:
    path "multiqc_report.html" into multiqc_output

    """
    multiqc $reports
    """
}

process mash_screen_genome {

    //errorStrategy 'ignore'
    publishDir params.out, mode: 'copy', overwrite: true

    input:
    tuple val(name), file(fastq) from fastq_files4

    output:
    path "*_pathogen_id.out" into mash_screen_genome_out

/*
    """
    cat ${fastq[0]} ${fastq[1]} > combined.fastq.gz
    mash screen -w ${mash_genome_db} combined.fastq.gz | sort -gr > ${name}_pathogen_id.out
    """ 
*/

    """
    cat ${fastq} | mash screen -w ${mash_genome_db} - | sort -gr > ${name}_pathogen_id.out
    """
}


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
}
*/

process kma_index {

    //errorStrategy 'ignore'
    //publishDir params.out, overwrite: true

    output:
    path "abr*" into kma_index_out

    """
    kma index -i ${abr_ref} -o abr
    """
}

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

process adapter_trimming {

    //errorStrategy 'ignore'
    //publishDir params.out, mode: 'copy', overwrite: true

    input:
    tuple val(name), file(fastq) from fastq_files3

    output:
    tuple val(name), file("*_clean{1,2}.fq.gz") into trimmed_fastq

    """
    bbduk.sh -Xmx1g in1=${fastq[0]} in2=${fastq[1]} out1=int1.fq.gz out2=int2.fq.gz ref=${phix} k=31 hdist=1 t=1
    bbduk.sh -Xmx1g in1=int1.fq.gz in2=int2.fq.gz out1=${name}_clean1.fq.gz out2=${name}_clean2.fq.gz ref=${adapters} ktrim=r k=23 mink=11 hdist=1 t=1 tpe tbo
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
    spades.py -1 ${fastq[0]} -2 ${fastq[1]} -o ${name} -m 16 -t ${params.thread}
    cp ${name}/scaffolds.fasta ${name}_scaffolds.fasta
    """
}

process assembly_size_filter {

    //errorStrategy 'ignore'
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

    //errorStrategy 'ignore'
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

    //errorStrategy 'ignore'
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

    errorStrategy 'ignore'
    publishDir params.out, mode: 'copy', overwrite: true

    input:
    tuple val(name), file(assembly) from assembly_filter_output3

    output:
    path("*") into busco_output

    """
    busco --auto-lineage-prok -f -m geno -o ${name}_busco -i ${assembly}
    """
}

