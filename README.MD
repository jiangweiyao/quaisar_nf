# QuAISAR Pipeline in Nextflow

### This project is still undercontruction.

### Dependencies:
1. [Nextflow](https://www.nextflow.io/)
2. [Singularity](https://sylabs.io/docs/) - if you are running this on a cluster, please install version 3.6 or later. 

### Usage: 
Test the full pipeline with prepackaged data:
```
nextflow run quaisar_nf/main.nf 
```

Use the configuration file to run on the SGE system in SCBS
```
nextflow run quaisar_nf/main.nf -c quaisar_nf/cdc_sge.config
```

Test the shorter KPIP pipeline
```
nextflow run quaisar_nf/kpip.nf 
```

Use the configuration file to run the shorter KPIP pipeline on the SGE system in SCBS
```
nextflow run quaisar_nf/kpip.nf -c quaisar_nf/cdc_sge.config
```

Put in your own parameters and get log tracing example:
```
nextflow run quaisar_nf/main.nf --in="quaisar_nf/fastq/*{1,2}.fastq.gz" --out="output_test" --abrdb="quaisar_nf/db/ARGannot/ARGannot_r3.fasta" --adapters="quaisar_nf/db/adapters/adapters.fa" --sizefilter=3000 -with-report report.html
```


