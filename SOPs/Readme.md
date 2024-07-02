# Setting up your Citizen Phage Environment

This pipeline is to be run in four stages:

## 1. Assemble your genome
This can be done using a command similar to 
```
snakemake -j 16 \
-s snakefiles/assembly.smk \
--use-conda \
output/reports/<phage_name>/assembly-report.json
```
This will create all the files needed to make a report of the assembly and output a report in JSON format

## 2. Manually select contigs for processing using Bandage
After lots of attempts to automate this process, nothing beats getting in there and manually checking and curating the assembly graphs to make sure you get the best virus. Single contigs should be saved in the folder `output/03_selected_contigs/` with the filename `<phage name>.fa` if there is only one of them or `<phage_name>_<number>.fa` if there is more than one (e.g. `CPL00001.fa` if there is only one contig or `CPL00001_1.fa` and `CPL00001_2.fa` if there are two). 

## 3. Check purity
It's always a good idea to check what is left over in the reads that map to neither your host or your viral contigs - this is a good way of finding contaminants, other phage remnants etc. This can be done using the command:
```
snakemake -j 16 \
-s snakefiles/purity-check.smk \
--use-conda \
output/reports/<phage_name>/purity-report.json
```

## 4. Run the annotation pipeline on each contig
This can be done using a command similar to:
```
snakemake -j 16 \
-s snakefiles/annotation.smk \
--use-conda \
output/reports/<phage_name>/<contig_name>/annotation-report.json

```
This will create all the files needed to make a report of the annotation and output a report in JSON format
