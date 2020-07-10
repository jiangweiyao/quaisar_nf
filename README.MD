# QuAISAR Pipeline in Nextflow

###This project is still undercontruction.

###Usage: 
Test pipeline with prepackaged data:
```
nextflow run quaisar_nf/main.nf
```

Put in your own parameters example:
```
nextflow/nextflow run quaisar_nf/main.nf --in="quaisar_nf/fastq/*{1,2}.fastq.gz" --out="output_test" --abrdb="quaisar_nf/db/ARGannot/ARGannot_r3.fasta" --adapters="quaisar_nf/db/adapter
s/adapters.fa" --sizefilter=3000 -resume
```
