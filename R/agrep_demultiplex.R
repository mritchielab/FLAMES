# work in progress
library("ShortRead")
library("Biostrings")
agrep_demultiplex_ <- function(fastq_file, bc.list, out.file = "out.fq.gz", structure = list(adaptor.p5 = "CTACAC",
  barcode = ".{16}", adaptor.read1nspacer = "CGCGTC.{28}GAGACAG", seq = ".*"),
  trim.after = "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC", strands = c("-", "+"), pattern.max.distance = 3,
  bc.max.distance = 2) {

  cat(fastq_file, "\n")
  strm <- FastqStreamer(fastq_file)
  nreads <- 0
  while (length(res <- yield(strm))) {
    nreads <- nreads + length(res)

    match.idx <- sapply(strands, function(x) {
      NULL
    })
    for (strand in sort(strands, decreasing = TRUE)) {
      if (strand == "-") {
        res <- reverseComplement(res)
      }

      pattern <- paste0(c("(", paste0(structure, collapse = ")("), ")"), collapse = "")
      m <- aregexec(pattern, sread(res), max.distance = pattern.max.distance)
      match.idx[[strand]] <- sapply(m, function(x) {
        x[1] != -1
      })
      res.matched <- res[match.idx[[strand]]]
      m <- m[match.idx[[strand]]]
      if (!length(m)) {
        next
      }

      if ("barcode" %in% names(structure)) {
        bc <- sapply(regmatches(sread(res.matched), m), function(x) {
          x[which("barcode" == names(structure)) + 1]
        })
          bc.exact <- bc %in% bc.list
        if (bc.max.distance >= 1 && length(bc[!bc.exact])) {
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
          bc.matched <- bc.exact
        }
        match.idx[[strand]][match.idx[[strand]]][!bc.matched] <- FALSE
        readseq <- DNAStringSet(sapply(regmatches(sread(res.matched), m),
          function(x) {
          x[which("seq" == names(structure)) + 1]
          }))[bc.matched]
        readid <- BStringSet(paste(bc, id(res.matched), sep = "#"))[bc.matched]
        readqual <- sapply(regmatches(quality(res.matched), m), function(x) {
          x[which("seq" == names(structure)) + 1]
        })[bc.matched]

      } else {
        readseq <- DNAStringSet(sapply(regmatches(sread(res.matched), m),
          function(x) {
          x[which("seq" == names(structure)) + 1]
          }))
        readid <- id(res.matched)
        readqual <- sapply(regmatches(quality(res.matched), m), function(x) {
          x[which("seq" == names(structure)) + 1]
        })
      }

      if (!missing(trim.after) && !is.null(trim.after)) {
        alig <- pairwiseAlignment(subject = trim.after, pattern = readseq,
          type = "overlap")
        readseq <- subseq(readseq, start = ifelse(score(alig) > 10, end(pattern(alig)),
          1))
        readqual <- SFastqQuality(subseq(readqual, start = ifelse(score(alig) >
          10, end(pattern(alig)), 1)))
      }

      demultiplexed_res <- new("ShortReadQ", sread = readseq, quality = readqual,
        id = readid)
      writeFastq(object = demultiplexed_res[width(demultiplexed_res) > 10],
        mode = "a", file = out.file)
    }
    if (length(strands) > 1) {
      cat(as.character(id(res[match.idx[["+"]] & match.idx[["-"]]])), sep = "\n")
      cat("\n")
    }
  }
  close(strm)

  stats <- data.frame(nreads = nreads)
  if (length(strands) > 1) {
    stats["n.mixed.strands"] <- sum(match.idx[["+"]] & match.idx[["-"]])
  }
  for (strand in strands) {
    stats[paste0("n.demultiplexed.", strand)] <- sum(match.idx[[strand]])
  }
  return(stats)
}

agrep_demultiplex <- function(fastq_file, bc.list, out.file = "out.fq.gz", structure = list(adaptor.p5 = "CTACAC",
  barcode = ".{16}", adaptor.read1nspacer = "CGCGTC.{28}GAGACAG", seq = ".*"),
  trim.after = "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC", strands = c("-", "+"), pattern.max.distance = 3,
  bc.max.distance = 2) {
  if (length(fastq_file) == 1) {
    return(agrep_demultiplex_(fastq_file, bc.list, out.file, structure, trim.after,
      strands, pattern.max.distance, bc.max.distance = 2))
  } else {
    bcs <- data.frame()
    for (i in fastq_file) {
      bcs <- rbind(bcs, agrep_demultiplex_(i, bc.list, out.file, structure,
        trim.after, strands, pattern.max.distance, bc.max.distance = 2))
    }
    return(bcs)
  }
}

barcodes.atac.tsv <- read.csv("../barcodes.atac.tsv", header = F)$V1
fq_files <- c("../pass/fastq_pass/PAM86628_pass_f70417de_14b6a3ab_1378.fastq.gz",
  "../pass/fastq_pass/PAM86628_pass_f70417de_14b6a3ab_221.fastq.gz", "../pass/fastq_pass/PAM86628_pass_f70417de_14b6a3ab_978.fastq.gz",
  "../pass/fastq_pass/PAM86628_pass_f70417de_14b6a3ab_311.fastq.gz", "../pass/fastq_pass/PAM86628_pass_f70417de_14b6a3ab_616.fastq.gz",
  "../pass/fastq_pass/PAM86628_pass_f70417de_14b6a3ab_23.fastq.gz", "../pass/fastq_pass/PAM86628_pass_f70417de_14b6a3ab_1117.fastq.gz",
  "../pass/fastq_pass/PAM86628_pass_f70417de_14b6a3ab_340.fastq.gz", "../pass/fastq_pass/PAM86628_pass_f70417de_14b6a3ab_498.fastq.gz",
  "../pass/fastq_pass/PAM86628_pass_f70417de_14b6a3ab_713.fastq.gz")

cat(format(Sys.time(), "%a %b %d %X"), "\n")
x <- agrep_demultiplex(fq_files, barcodes.atac.tsv, bc.max.distance = 1, pattern.max.distance = 3)
cat(format(Sys.time(), "%a %b %d %X"), "\n")

saveRDS(x, "stats.rds")


# structure <- list(`-` = list(seq = '.*', adaptor.read1nspacer =
# 'CTGTCTC.{28}GACGCG', barcode = '.{16}', adaptor.p5 = 'GTGTAG'), `+` =
# list(adaptor.p5 = 'CTACAC', barcode = '.{16}', adaptor.read1nspacer =
# 'CGCGTC.{28}GAGACAG', seq = '.*'))

#'((\\[.+\\])|([A-Za-z\\.]))(\\{[0-9,]+\\})'
# complement_regex <- function(str) { new_str = '' for (i in strsplit(str,
# split = '')) { new_str <- paste0(c(new_str, switch(i, A = 'T', C = 'G', T =
# 'A', G = 'C', i))) return(new_str) } } reverse_pattern <- '(' for
# (regex_group in rev(structure)) { }
