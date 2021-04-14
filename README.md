# findmotifs
An R script for finding and scoring motifs in small genomes using position-frequency tables.  

### Inputs

1. fasta file (`-f`,`--fasta`)
2. position-frequency table (`-t`,`--freq_table`)
3. minimum score for matches (`-c`,`--cutoff`; default is 90)
4. output directory (`-o`,`--out_dir`; default is `getwd()`)

The position-frequency table must be a concatenated list of comma-delimited position-frequency matrices, with the name of each motif preceding its corresponding matrix. For an example, see the structure of `PFMs.txt`: 

    motif_A
    A,58,46,41,0,0,28,34,0,23,34,0,64
    C,9,48,9,11,95,0,0,92,28,51,60,2
    G,33,0,24,0,5,72,47,8,17,15,0,0
    T,0,6,26,89,0,0,19,0,32,0,40,34
    motif_B
    A,0,48,8,0,0,100,37,28,23,0,29,1,0
    C,100,0,35,24,0,0,0,9,34,100,71,0,0
    G,0,52,41,16,100,0,63,36,35,0,0,0,100
    T,0,0,16,60,0,0,0,27,8,0,0,99,0
    motif_C
    A,0,0,93,0,31,25,0,0,0,29
    C,0,71,0,0,0,52,0,0,100,59
    G,7,0,0,47,69,23,90,0,0,0
    T,93,29,7,53,0,0,10,100,0,12

### Outputs

This script returns a single tab-delimited file with the following fields:

`motif_name`  
- the name of the motif as it appears in the position-frequency table  

`motif_width`  
- the number of nucleotides in the motif  

`fasta`  
- the name of the input fasta  

`fasta_length`  
- the number of nucleotides in the input fasta  

`strand`  
- the strand the match is found on (`+` or `-`)  

`match`  
- the nucleotide sequence of the match  

`score`  
- the score of the match  

`concat_start`  
- the position of the first nucleotide of the match in the `+` concatenated fasta file  

`concat_end`  
- the position of the last nucleotide of the match in the `+` concatenated fasta file. For matches on the `+` strand, `concat_start` < `concat_end`. For matches on the `-` strand, `concat_start` > `concat_end`.  

`contig_index`  
- the contig number on which the match is located  

`contig_length`  
- the number of nucleotides in the contig containing the match  

`contig_start`  
- the position of the first nucleotide in the match relative to the `+` start of the contig  

`contig_end`  
- the position of the last nucleotide in the match relative to the `+` start of the contig. For matches on the `+` strand, `contig_start` < `contig_end`. For matches on the `-` strand, `contig_start` > `contig_end`.  

`contig_header`  
- the fasta header of the contig

### Dependencies

1. R, version >= 4.0.0
2. [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html) (tested version: 2.56.0)
3. [optparse](https://cran.r-project.org/web/packages/optparse/index.html) (tested version: 1.6.6)
4. [stringr](https://cran.r-project.org/web/packages/stringr/index.html) (tested version: 1.4.0)

### Example

To run `findmotifs.R` from the command line, ensure R is in your `PATH` and call the script with `Rscript`. 

    export PATH=/programs/R-4.0.0/bin:$PATH  
    Rscript findmotifs.R -f EcoliC.fa -t PFMs.txt -c 90 -o .

Use `-h` or `--help` to print the help message.

    Rscript findmotifs.R -h
    Options:
        -f CHARACTER, --fasta=CHARACTER
                Required: input directory and filename of fasta sequence

        -t CHARACTER, --freq_table=CHARACTER
                Required: input directory and filename of comma-delimited position frequency table

        -c CHARACTER, --cutoff=CHARACTER
                Input cutoff as integer <= 100 (default: 90)

        -o CHARACTER, --out_dir=CHARACTER
                Output directory without trailing "/" (default: getwd())

        -h, --help
                Show this help message and exit

### Details

This script is a wrapper for the `matchPWM()` function from the `Biostrings` package. (see [RDocumentation](https://www.rdocumentation.org/packages/Biostrings/versions/2.40.2/topics/matchPWM))

To save time, the number of `matchPWM()` calls are reduced to one per motif by first concatenating the input fasta into a single sequence with contigs demarcated by non-IUPAC characters. 

Future implementations should further increase speed by processing match tables as data frames and not hit-by-hit. 

