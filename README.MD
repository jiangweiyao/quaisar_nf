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

### Resuming the workflow
If your workflow is interrupted, you can resume it by running the same command with -resume after it to reuse the already successfully executed steps. For example:
```
nextflow run quaisar_nf/kpip.nf -resume
```


### Common Input Patterns:
Files will be searched against the glob expression for the input argument (i.e. --in="quaisar_nf/fastq/*{1,2}.fastq.gz").
Common patterns include:
1. --in="<path>/*{1,2}.fastq.gz" for paired SRA files
2. --in="<path>/*{R1_001,R2_001}.fastq.gz" for paired Illumina files

You can find all files in subfolders of a directory using the following expression:
--in="<path>/**/*{R1_001,R2_001}.fastq.gz"

You can determine which files match the glob pattern using the include filepair_finder.nf workflow before you run the workflow on your actual files. You can try running this against the included fastq files using the following command:
```
nextflow run quaisar_nf/filepair_finder.nf --in='./quaisar_nf/fastq/*{1,2}.fastq.gz'
```

### Common problems
This workflow will dynamically download the docker images for the processes and the large reference files at first running. Network interruptions or firewalls can block the download. We recommend running the workflow with the include file the first time you run this workflow to download everything. If a docker image pull is unsuccessful, try resuming the run again to see if the docker image can be downloaded successfully. 
