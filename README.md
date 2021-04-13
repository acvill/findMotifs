# findmotifs
An R script for finding and scoring motifs in small genomes using position-frequency tables

### Inputs

1. A fasta file (`-f`, `--fasta`)
2. 
3. 

### Outputs


### Dependencies

1. R, version >= 4.0.0
2. [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html) (tested version: 2.56.0)
3. [optparse](https://cran.r-project.org/web/packages/optparse/index.html) (tested version: 1.6.6)
4. [stringr](https://cran.r-project.org/web/packages/stringr/index.html) (tested version: 1.4.0)

### Example

    export PATH=/programs/R-4.0.0/bin:$PATH  
    Rscript findmotifs.R -f EcoliC.fa -t PFMs.csv -c 90% -o .

Use `-h` or `--help` to print the help message

    Rscript findmotifs.R -h
    Options:
        -f CHARACTER, --fasta=CHARACTER
                Required: input directory and filename of fasta sequence

        -t CHARACTER, --freq_table=CHARACTER
                Required: input directory and filename of comma-delimited position frequency table

        -c CHARACTER, --cutoff=CHARACTER
                Input cutoff as 'integer%' (default: 90%)

        -o CHARACTER, --out_dir=CHARACTER
                Output directory without trailing "/" (default: getwd())

        -h, --help
                Show this help message and exit
