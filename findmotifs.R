#!/usr/bin/env Rscript

# last edited 14 April 2021, Albert Vill, acv46@cornell.edu

# Functionality
## concatenates input fastas into a single DNA string with non-DNA characters at contig boundaries
## the concat_fasta is then searched using each position weight matrix and matches are parsed to get contig-wise start and end positions

# load libraries
suppressMessages(library("Biostrings"))
suppressMessages(library("optparse"))
suppressMessages(library("stringr"))

# make_pfm code to generate mock frequency matrices
# rand_vect function from https://stackoverflow.com/a/24846002/7976890

# rand_vect <- function(N, M, sd = 1, pos.only = TRUE) {
#   vec <- rnorm(N, M/N, sd)
#   if (abs(sum(vec)) < 0.01) vec <- vec + 1
#   vec <- round(vec / sum(vec) * M)
#   deviation <- M - sum(vec)
#   for (. in seq_len(abs(deviation))) {
#     vec[i] <- vec[i <- sample(N, 1)] + sign(deviation)
#   }
#   if (pos.only) while (any(vec < 0)) {
#     negs <- vec < 0
#     pos  <- vec > 0
#     vec[negs][i] <- vec[negs][i <- sample(sum(negs), 1)] + 1
#     vec[pos][i]  <- vec[pos ][i <- sample(sum(pos ), 1)] - 1
#   }
#   vec
# }
# 
# make_pfm <- function(x, M = 100, sd = 20) {
#   for (i in 1:x) {
#     if (i == 1) {
#       mat <- matrix(rand_vect(4, M = M, sd = sd))
#     }
#     else {
#       mat <- cbind(mat, matrix(rand_vect(4, M = M, sd = sd)))
#     }
#   }
#   row.names(mat) <- c("A","C","G","T")  
#   mat
# }

# handle input as command-line options with optparse
option_list=list(
  make_option(c("-f", "--fasta"),
              action = "store",
              type = "character",
              default = NULL,
              help = "REQUIRED: input directory and filename of fasta sequence",
              metavar = "character"),
  make_option(c("-t", "--freq_table"),
              action = "store",
              type = "character",
              default = NULL,
              help = "REQUIRED: input directory and filename of comma-delimited position frequency table",
              metavar = "character"),
  make_option(c("-c", "--cutoff"),
              action = "store",
              type = "integer",
              default = "90",
              help = "Input cutoff as integer <= 100 (default: 90)",
              metavar = "character"),
  make_option(c("-o", "--out_dir"),
              action = "store",
              type = "character",
              default = getwd(),
              help = "Output directory without trailing \"/\" (default: getwd())",
              metavar = "character"))
opt_parser=OptionParser(option_list = option_list)
opt=parse_args(opt_parser)

message(paste("## reading position frequency table -->", opt$freq_table, sep = " "))
message(paste("## reading fasta file -->", opt$fasta, sep = " "))
message(paste("## score cutoff set at ", opt$cutoff, "%", sep = ""))

# create temp folder for intermediate files
if (dir.exists(paths = paste(opt$out_dir)) == T){
  rand <- as.numeric(paste(sample(9,10,replace = T), collapse = ""))
  tpath <- paste(opt$out_dir, "/temp", rand, sep = "")
  dir.create(path = tpath, showWarnings = TRUE, recursive = FALSE, mode = "0777")
} else if (dir.exists(paths = paste(opt$out_dir)) == F){
  message("The supplied --out_dir does not exist. Aborting program.")
  stop()
}

message(paste("## output directory set to", opt$out_dir, sep = " "))

# read in position frequency table object from which position weight matrices will be generated
numcol <- max(count.fields(file = paste(opt$freq_table),
                           sep = ","))
linear <- read.table(file = paste(opt$freq_table),
                     header=F,
                     sep=",",
                     fill=T,
                     col.names=paste0('V', seq_len(numcol)))

# read in fasta as DNAStringSet object and get summary statistics
fasta <- readDNAStringSet(filepath = opt$fasta,
                          format = "fasta")
numseq <- length(fasta)
fname <- basename(opt$fasta)
message(paste(fname, "contains", numseq, "contigs ... concatenating", sep = " "))


# concatenate contigs into single nucleotide string
# contigs separated by periods, which is a valid character in Biostrings::DNA_ALPHABET
concat_fasta <- c()
for (i in 2:numseq){
  
  if (i == 2){
    
    concat_fasta <- paste(as.character(fasta[[i-1]]),
                          c(".........."),
                          as.character(fasta[[i]]),
                          sep = '',
                          collapse = '')
    
  }else{
    
    concat_fasta <- paste(concat_fasta,
                          c(".........."),
                          as.character(fasta[[i]]),
                          sep = '',
                          collapse = '')
    
  }
  
}

# do the same as above, but for the reverse complement
concat_fasta_rev <- c()
for (i in 1:(numseq - 1)){
  
  if (i == 1){
    
    concat_fasta_rev <- paste(as.character(reverseComplement(fasta[[numseq]])),
                          c(".........."),
                          as.character(reverseComplement(fasta[[numseq - i]])),
                          sep = '',
                          collapse = '')
    
  }else{
    
    concat_fasta_rev <- paste(concat_fasta_rev,
                          c(".........."),
                          as.character(reverseComplement(fasta[[numseq - i]])),
                          sep = '',
                          collapse = '')
    
  }
  
}

count1 <- as.integer(nchar(concat_fasta) - ((numseq - 1)*10))
count2 <- sum(fasta.seqlengths(filepath = opt$fasta))
count3 <- as.integer(nchar(concat_fasta_rev) - ((numseq - 1)*10))
dotcount <- sum(str_count(fasta, "\\."))

# Check if concatentated fasta has the same ATGC character length as original fasta
if (count1 != count2){
  
  message(paste("ERROR loading", fname, "-- length of concatenated fasta does not agree with Biostrings::fasta.seqlengths(). Program aborted.",
                sep =  ' '))
  stop()
  
}

# Check if reverse concatenated fasta has the same ATGC character length as original fasta
if (count1 != count3){
  
  message(paste("ERROR loading", fname, "-- different lengths for forward and reverse fastas. Check ?Biostrings::reverseComplement() for help. Program aborted.",
                sep =  ' '))
  stop()
  
}

# Check for the absence of period characters in the sequence content of the input fasta
if (dotcount != 0){
  
  message(paste("ERROR loading", fname, "-- fasta sequence contains", dotcount, "periods (\".\"), which are a special character used by this script. Program aborted.",
                sep =  ' '))
  stop()
  
}

# read in motif list as vector
# assumes every fifth line is a motif identifier
motif_list <- linear[seq(1, nrow(linear), 5), 1]

# generate column names object for results table
rescols <- c("motif_name",
             "motif_width",
             "fasta",
             "fasta_length",
             "strand",
             "match",
             "score",
             "concat_start",
             "concat_end",
             "contig_index",
             "contig_length",
             "contig_start",
             "contig_end",
             "contig_header")

# calculate GC content of fasta to adjust priors of PWM function
gcc <- round(sum(letterFrequency(x = fasta, letters = c("G","C"))) / count1, 4)
prior <- c(A=(1-gcc)/2, C=gcc/2, G=gcc/2, T=(1-gcc)/2)

for (motif in 1:length(motif_list)){
 
  # generate rbp object to store motif name and row index of motif in "linear" object
  rbp <- as.character(append(motif_list[motif], as.numeric(which(linear$V1 == motif_list[motif]))))

  # for a given RBP, generate pwm from input position-weight data
  pfm <- linear[(as.numeric(rbp[2])+1):(as.numeric(rbp[2])+4),]
  pfm <- pfm[,colSums(is.na(pfm)) == 0]
  row.names(pfm) <- pfm[,1]
  pfm <- data.matrix(pfm[,-1])
  pwm <- PWM(x = pfm,
             prior.params = prior,
             type = c("log2probratio"))
  
  message(paste("Position weight matrix generated for motif", rbp[1], sep = " "))
  
  # find matches on forward strand
    
  fwdmatch <- suppressWarnings(matchPWM(pwm = pwm,
                                        subject = concat_fasta,
                                        min.score = paste(opt$cutoff,"%"),
                                        with.score = T))
  
  numfwd <- length(fwdmatch)
  message(paste(rbp[1], "--", numfwd, "matches found on fwd strand of", fname, sep = " "))
  
  # find matches for reverse strand
  
  revmatch <- suppressWarnings(matchPWM(pwm = pwm,
                                        subject = concat_fasta_rev,
                                        min.score = paste(opt$cutoff,"%"),
                                        with.score = T))
  
  numrev <- length(revmatch)
  message(paste(rbp[1], "--", numrev, "matches found on rev strand of", fname, sep = " "))
  
  message(paste("Processing matches for", rbp[1], sep = " "))
  
  # initiate progress bar
  pbar = txtProgressBar(min = 0,
                        max = numfwd+numrev,
                        initial = 0,
                        style = 3,
                        file = "") 

  # If matches are found on plus / fwd strand, write results to object
  if (numfwd > 0){
  
    fwdres <- data.frame()
    for (mch in 1:numfwd){

      # get number of period characters between the start of concat_fasta and the start of the match
      dot_index <- str_count(substr(concat_fasta, 1, start(fwdmatch[mch,1])), "\\.")
      # get the contig context of the match, based on count of period spacers
      contig_index <- (dot_index / 10) + 1
      # for plus strand concat_start < concat_stop -- see note for minus strand annotation below
      concat_start <- start(fwdmatch[mch,1]) - dot_index
      motif_width <- width(fwdmatch[mch,1])
      concat_end <- concat_start + motif_width - 1
      # get the contig-wise start and end positions of the match, match width, contig length, and contig header
      contig_length <- width(fasta[contig_index])
      contig_header <- names(fasta[contig_index])
      contig_start <- concat_start - sum(width(fasta[0:(contig_index-1)]))
      contig_end <- contig_start + motif_width - 1
      
      fwdpart <- c(motif_name = rbp[1],
                   motif_width = motif_width,
                   fasta = fname,
                   fasta_length = count1,
                   strand = "+",
                   match = data.frame(fwdmatch)[,1][mch],
                   score = round(mcols(fwdmatch)$score[mch],5),
                   concat_start = concat_start,
                   concat_end = concat_end,
                   contig_index = contig_index,
                   contig_length = contig_length,
                   contig_start = contig_start,
                   contig_end = contig_end,
                   contig_header = contig_header)
    
      # append match line to object representing all matches for motif
      fwdres <- rbind(fwdres, fwdpart)
      
      # iterate progress bar
      setTxtProgressBar(pb = pbar,
                        value = mch)
      
    }
  
    colnames(fwdres) <- rescols  

  }
  
  # If matches are found on minus / rev strand, write results to object
  
  if (numrev > 0){
  
    revres <- data.frame()
    for (mch in 1:numrev){
    
      # get number of period characters between the start of concat_fasta_rev and the start of the match
      dot_index <- str_count(substr(concat_fasta_rev, 1, start(revmatch[mch,1])), "\\.")
      # get the contig context of the match, based on count of period spacers
      contig_index <- numseq - (dot_index / 10)
      # for minus strand, start and end positions given relative to plus strand, such that concat_start > concat_stop
      concat_end <- count2 - (end(revmatch[mch,1]) - ((numseq - contig_index) * 10))
      motif_width <- width(revmatch[mch,1])
      concat_start <- concat_end + motif_width - 1
      # get the contig-wise start and end positions of the match, match width, contig length, and contig header
      contig_length <- width(fasta[contig_index])
      contig_header <- names(fasta[contig_index])
      contig_end <- concat_end - sum(width(fasta[0:(contig_index-1)]))
      contig_start <- contig_end + motif_width - 1
      
      revpart <- c(motif_name = rbp[1],
                   motif_width = motif_width,
                   fasta = fname,
                   fasta_length = count1,
                   strand = "-",
                   match = data.frame(revmatch)[,1][mch],
                   score = round(mcols(revmatch)$score[mch],5),
                   concat_start = concat_start,
                   concat_end = concat_end,
                   contig_index = contig_index,
                   contig_length = contig_length,
                   contig_start = contig_start,
                   contig_end = contig_end,
                   contig_header = contig_header)
      
      # append match line to object representing all matches for motif
      revres <- rbind(revres, revpart)
      
      # iterate progress bar
      setTxtProgressBar(pb = pbar,
                        value = numfwd + mch)
    
    }
  
    colnames(revres) <- rescols
  
  }

  # write object containing forward and reverse hits to temporary file
  # combining objects with rbind breaks for objects with zero lines
  # check that both fwd and rev hits are nonzero before combining
  if (numfwd > 0 && numrev > 0){
    
    message("")
    message(paste("Appending forward and reverse matches to" , rbp[1], "to results table for", fname, sep = " "))
    seqres <- rbind(fwdres, revres)
    
  } else if (numfwd > 0 && numrev == 0) {
    
    message("")
    message(paste("No hits on minus strand for", rbp[1], ", appending forward hits to results table for ", fname, sep = ""))
    seqres <- fwdres
    
  } else if (numfwd == 0 && numrev > 0) {
    
    message("")
    message(paste("No hits on plus strand for", rbp[1], ", appending forward hits to results table for ", fname, sep = ""))
    seqres <- revres
    
  } else if (numfwd == 0 && numrev == 0) {  
    
    message("")
    message(paste("No hits on either strand for ", rbp[1], " in ", fname, ". Try reducing cutoff.", sep = ""))
    stop()
    
  }
  
  write.table(seqres,
              file = paste(tpath,
                           "/",
                           opt$fasta,
                           "_",
                           motif_list[motif],
                           ".temp",
                           sep = ""),
              sep = ",")
  
}

# combine temporary files into single file containing hits for all input motifs

setwd(tpath)
mergedhits <- do.call(rbind,
                      lapply(list.files(path = tpath), read.csv))

setwd(opt$out_dir)
write.table(x = mergedhits, 
            row.names=F, 
            sep="\t", 
            quote = F,
            file = paste(opt$out_dir, paste(opt$fasta, "_findmotifs.txt", sep = ""), sep = "/"))

message(paste("Removing temporary files"))
unlink(tpath, recursive = T)    
