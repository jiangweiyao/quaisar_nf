# QuAISAR Pipeline in Nextflow

###This project is still undercontruction.

###Usage: 
Test pipeline with prepackaged data:
```
nextflow run quaisar_nf/main.nf -c quaisar_nf/nextflow.config
```

Put in your own parameters and get log tracing example:
```
nextflow run quaisar_nf/main.nf -c quaisar_nf/nextflow.config --in="quaisar_nf/fastq/*{1,2}.fastq.gz" --out="output_test" --abrdb="quaisar_nf/db/ARGannot/ARGannot_r3.fasta" --adapters="quaisar_nf/db/adapters/adapters.fa" --sizefilter=3000 -with-report report.html
```
