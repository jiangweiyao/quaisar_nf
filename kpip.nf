#!/usr/bin/env nextflow


Channel.fromFilePairs( params.in ).into { fastq_files; fastq_files2; fastq_files3; fastq_files4; fastq_files5; fastq_files6; fastq_files7  }
abr_ref = file(params.abrdb)
plasmid_db = file(params.plasmid_db)
adapters = file(params.adapters)
phix = file(params.phix)
params.thread = 1
mash_genome_db = file(params.genome_db)
kraken_db = file(params.kraken_db)

busco_config = file("$baseDir/db/busco_config.ini")
mash_parser = file("$baseDir/bin/mash_screen_parser.py")
kma_parser = file("$baseDir/bin/kma_screen_parser.py")

println """\
         K P I P     NEXTFLOW      PIPELINE   
         =====================================================
         input reads (--in)                  : ${params.in}
         outdir (--out)                      : ${params.out}
         Antibiotic Resistance DB (--abrdb)  : ${params.abrdb}
         Mash Genome Reference (--genome_db) : ${params.genome_db}
         Plasmid Database (--plasmid_db)     : ${params.plasmid_db}

         """
         .stripIndent()


// Check if Mash reference file already exists. If not, download it.
mash_genome_file = file("$baseDir/db/refseq.genomes.k21s1000.msh")
if(!mash_genome_file.exists()){
    println("Mash genome reference missing. Downloading...")
    mash_genome_file = file('https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh')
    mash_genome_file.copyTo("$baseDir/db/refseq.genomes.k21s1000.msh")
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
    tuple val(name), path("*_pathogen_id.out") into mash_screen_genome_out

    """
    cat ${fastq} | mash screen -w ${mash_genome_db} - | sort -gr > ${name}_pathogen_id.out
    """
}

process tabulate_mash_genome {

    //errorStrategy 'ignore'
    publishDir params.out, mode: 'copy', overwrite: true

    input:
    tuple val(name), file(table) from mash_screen_genome_out

    output:
    path "*_pathogen_id.out1" into tabulate_mash_genome_out

    """
    python3 ${mash_parser} -i ${table} -o ${name}_pathogen_id.out1
    """
}


process kma_index_abr {

    //errorStrategy 'ignore'
    //publishDir params.out, overwrite: true

    output:
    path "abr*" into kma_index_abr_out

    """
    kma index -i ${abr_ref} -o abr
    """
}

process kma_map_abr {

    //errorStrategy 'ignore'
    publishDir params.out, mode: 'copy', overwrite: true

    input:
    path index from kma_index_abr_out
    tuple val(name), file(fastq) from fastq_files6

    output:
    tuple val(name), file("*_abr*") into kma_abr_out

    //kma -i ${fastq[0]} -ipe ${fastq[1]} -o ${name}_abr.out -t_db abr -1t1


    """
    kma -ipe ${fastq} -o ${name}_abr -t_db abr -1t1
    """
}

process tabulate_kma_abr {

    //errorStrategy 'ignore'
    publishDir params.out, mode: 'copy', overwrite: true

    input:
    tuple val(name), file(table) from kma_abr_out

    output:
    path "*_abr.out1" into tabulate_kma_abr_out

    """
    python3 ${kma_parser} -i ${table[3]} -o ${name}_abr.out1
    """
}

process kma_index_plasmid {

    //errorStrategy 'ignore'
    //publishDir params.out, overwrite: true

    output:
    path "plasmid*" into kma_index_plasmid_out

    """
    kma index -i ${plasmid_db} -o plasmid
    """
}

process kma_map_plasmid {

    //errorStrategy 'ignore'
    publishDir params.out, mode: 'copy', overwrite: true

    input:
    path index from kma_index_plasmid_out
    tuple val(name), file(fastq) from fastq_files2

    output:
    tuple val(name), file("*_plasmid*") into kma_plasmid_out


    """
    kma -ipe ${fastq} -o ${name}_plasmid -t_db plasmid -1t1
    """
}

process tabulate_kma_plasmid {

    //errorStrategy 'ignore'
    publishDir params.out, mode: 'copy', overwrite: true

    input:
    tuple val(name), file(table) from kma_plasmid_out

    output:
    path "*_plasmid.out1" into tabulate_kma_plasmid_out

    """
    python3 ${kma_parser} -i ${table[3]} -o ${name}_plasmid.out1
    """
}


