# work in progress
library("ShortRead")
library("Biostrings")
agrep_demultiplex_ <- function(fastq_file, bc.list, out.file = "out.fq.gz", pattern = "(.*)(CTGTCTC.{28}GACGCG)(.{16})(GTGTAG)",
  pattern.max.distance = 3, bc.max.distance = 2, trim.after = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG") {
  bcs <- c()
  strm <- FastqStreamer(fastq_file)
  while (nreads <- length(res <- yield(strm))) {
    m <- aregexec(pattern, sread(res), max.distance = pattern.max.distance)

    res <- res[sapply(m, function(x) {
      x[1] != -1
    })]
    m <- m[sapply(m, function(x) {
      x[1] != -1
    })]

    bc <- sapply(regmatches(sread(res), m), function(x) {
      x[4]
    })

    if (bc.max.distance >= 1) {
      bc.exact <- bc %in% bc.list
      bc.dists <- sapply(bc[!bc.exact], function(x) {
        adist(x, bc.list)
      })
      bc.fuzzy <- apply(bc.dists, 2, function(x) {
        ifelse(min(x) < bc.max.distance && sum(x == min(x)) == 1, which.min(x),
          -1)
      })
      bc[!bc.exact][bc.fuzzy != -1] <- bc.list[bc.fuzzy[bc.fuzzy != -1]]

      bc.matched <- bc.exact
      bc.matched[!bc.exact][bc.fuzzy != -1] <- TRUE
    } else {
      bc.matched <- bc %in% bc.list
    }

    readseq <- DNAStringSet(sapply(regmatches(sread(res), m), function(x) {
      x[2]
    }))[bc.matched]
    readid <- BStringSet(paste(bc, id(res), sep = "#"))[bc.matched]

    alig <- pairwiseAlignment(subject = trim.after, pattern = readseq, type = "overlap")
    readseq <- subseq(readseq, start = ifelse(score(alig) > 10, end(pattern(alig)),
      1))
    readqual <- sapply(regmatches(quality(res), m), function(x) {
      x[2]
    })[bc.matched]

    readqual <- SFastqQuality(subseq(readqual, start = ifelse(score(alig) > 10,
      end(pattern(alig)), 1)))

    demultiplexed_res <- new("ShortReadQ", sread = readseq, quality = readqual,
      id = readid)
    writeFastq(object = demultiplexed_res[width(demultiplexed_res) > 10], mode = "a",
      file = out.file)

    length(bc) <- nreads
    bcs <- append(bcs, bc)
  }

  close(strm)
  return(bcs)
}

agrep_demultiplex <- function(fastq_file, bc.list, out.file = "out.fq.gz", pattern = "(.*)(CTCTGTC.{28}GACGCG)(.{16})(GTGTAG)",
  pattern.max.distance = 3, bc.max.distance = 2, trim.after = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG") {
  if (length(fastq_file) == 1) {
    return(agrep_demultiplex_(fastq_file, bc.list, out.file, pattern, pattern.max.distance,
      bc.max.distance = 2, trim.after))
  } else {
    bcs <- c()
    for (i in fastq_file) {
      bcs <- append(bcs, agrep_demultiplex_(i, bc.list, out.file, pattern,
        pattern.max.distance, bc.max.distance, trim.after))
    }
    return(bcs)
  }
}

barcodes.atac.tsv <- as.character(reverseComplement(DNAStringSet(read.csv("../barcodes.atac.tsv",
  header = F)$V1)))
fq_files <- c("../pass/fastq_pass/PAM86628_pass_f70417de_14b6a3ab_1378.fastq.gz",
  "../pass/fastq_pass/PAM86628_pass_f70417de_14b6a3ab_221.fastq.gz", "../pass/fastq_pass/PAM86628_pass_f70417de_14b6a3ab_978.fastq.gz",
  "../pass/fastq_pass/PAM86628_pass_f70417de_14b6a3ab_311.fastq.gz", "../pass/fastq_pass/PAM86628_pass_f70417de_14b6a3ab_616.fastq.gz",
  "../pass/fastq_pass/PAM86628_pass_f70417de_14b6a3ab_23.fastq.gz", "../pass/fastq_pass/PAM86628_pass_f70417de_14b6a3ab_1117.fastq.gz",
  "../pass/fastq_pass/PAM86628_pass_f70417de_14b6a3ab_340.fastq.gz", "../pass/fastq_pass/PAM86628_pass_f70417de_14b6a3ab_498.fastq.gz",
  "../pass/fastq_pass/PAM86628_pass_f70417de_14b6a3ab_713.fastq.gz")

cat(format(Sys.time(), "%a %b %d %X"), "\n")
x <- agrep_demultiplex(fq_files, barcodes.atac.tsv, bc.max.distance = 1, pattern.max.distance = 3)
cat(format(Sys.time(), "%a %b %d %X"), "\n")

#saveRDS(x, "bcs.rds")


"(.*)(CTCTGTC.{28}GACGCG)(.{16})(GTGTAG)"
structure <- list(seq = ".*", adaptor1 = "CTCTGTC.{28}GACGCG", barcode = ".{16}",
  adaptor2 = "GTGTAG")
forward_pattern <- paste0(c("(", paste0(structure, collapse = ")("), ")"), collapse = "")
"((\\[.+\\])|([A-Za-z\\.]))(\\{[0-9,]+\\})"

complement_regex <- function(str) {
  new_str = ""
  for (i in strsplit(str, split = "")) {
    new_str <- paste0(c(new_str, switch(i, A = "T", C = "G", T = "A", G = "C",
      i)))
    return(new_str)
  }
}

reverse_pattern <- '('
for (regex_group in rev(structure)) {
  reverse_pattern <- paste0()
}
