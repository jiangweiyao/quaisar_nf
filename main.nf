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
         Q U A I S A R - H     NEXTFLOW      PIPELINE   
         =====================================================
         input reads (--in)                  : ${params.in}
         outdir (--out)                      : ${params.out}
         Antibiotic Resistance DB (--abrdb)  : ${params.abrdb}
         adapters (--adapters)               : ${params.adapters}
         Mash Genome Reference (--genome_db) : ${params.genome_db}
         size cutoff (--sizefilter)          : ${params.sizefilter}
         Plasmid Database (--plasmid_db)     : ${params.plasmid_db}
         Kraken Database (--kraken_db)       : ${params.kraken_db}

         """
         .stripIndent()


// Check if Mash reference file already exists. If not, download it.
mash_genome_file = file("$baseDir/db/refseq.genomes.k21s1000.msh")
if(!mash_genome_file.exists()){
    println("Mash genome reference missing. Downloading...")
    mash_genome_file = file('https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh')
    mash_genome_file.copyTo("$baseDir/db/refseq.genomes.k21s1000.msh")
}

// Check if kraken2 library already exists. If not, download it.
kraken_hash_file = file("$baseDir/db/minikraken2_v2_8GB_201904_UPDATE/hash.k2d")
if(!kraken_hash_file.exists()){
    println("Kraken library missing. Downloading...")
    kraken_file = file('ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/old/minikraken2_v2_8GB_201904.tgz')
    //kraken_file = file('ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/16S_Greengenes13.5_20200326.tgz')
    kraken_file.copyTo("${baseDir}/db/minikraken2_v2_8GB_201904.tgz")
    //println "wget ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/old/minikraken2_v2_8GB_201904.tgz -P ${baseDir}/db/".execute().text
    //println "echo untarring tar -zxf ${baseDir}/db/minikraken2_v2_8GB_201904.tgz".execute().text
    println "tar -zxf ${baseDir}/db/minikraken2_v2_8GB_201904.tgz --directory ${baseDir}/db/".execute().text
    println "rm -f ${baseDir}/db/minikraken2_v2_8GB_201904.tgz".execute().text
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

process kraken_fastq {

    //errorStrategy 'ignore'
    publishDir params.out, mode: 'copy', overwrite: true

    input:
    tuple val(name), file(fastq) from fastq_files7

    output:
    tuple val(name), file("*.{summary,output}") into kraken_fastq_out

    """
    kraken2 --db ${kraken_db} --paired ${fastq} --memory-mapping --report ${name}_reads.summary --output ${name}_reads.output
    """
}

process krona_fastq {

    //errorStrategy 'ignore'
    publishDir params.out, mode: 'copy', overwrite: true

    input:
    tuple val(name), file(kraken_result) from kraken_fastq_out

    output:
    path "*.html" into krona_fastq_output

    """
    ktImportTaxonomy -q 2 -t 3 ${kraken_result[0]} -o ${name}_read.html
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
    tuple val(name), path("*_scaffolds_filtered.fasta") into assembly_filter_output, assembly_filter_output2, assembly_filter_output3, assembly_filter_output4, assembly_filter_output5

    """
    reformat.sh in=${assembly} out=${assembly.simpleName}_filtered.fasta minlength=${params.sizefilter}
    """
}

process kraken_assembly {

    //errorStrategy 'ignore'
    publishDir params.out, mode: 'copy', overwrite: true

    input:
    tuple val(name), file(assembly) from assembly_filter_output5

    output:
    tuple val(name), file("*.{summary,output}") into kraken_assembly_out

    """
    kraken2 --db ${kraken_db} ${assembly} --memory-mapping --report ${name}_assembly.summary --output ${name}_assembly.output
    """
}

process krona_assembly {

    //errorStrategy 'ignore'
    publishDir params.out, mode: 'copy', overwrite: true

    input:
    tuple val(name), file(kraken_result) from kraken_assembly_out

    output:
    path "*.html" into krona_assembly_output

    """
    ktImportTaxonomy -q 2 -t 3 ${kraken_result[0]} -o ${name}_assembly.html
    """
}


process mlst {

    //errorStrategy 'ignore'
    publishDir params.out, mode: 'copy', overwrite: true

    input:
    tuple val(name), file(assembly) from assembly_filter_output4

    output:
    tuple val(name), path("*.mlst") into mlst_output

    """
    mlst ${assembly} > ${name}.mlst
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
    busco --auto-lineage-prok -f -m geno -o ${name}_busco -i ${assembly} --config ${busco_config}
    """
}

