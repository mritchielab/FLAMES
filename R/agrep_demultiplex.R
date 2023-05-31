# work in progress
agrep_demultiplex_ <- function(fastq_file, bc.list, out.file, structure, trim.after,
  strands, pattern.max.distance, bc.max.distance) {
  cat(fastq_file, "\n")
  stat.names <- c("matched", "found.2N", "bc", "bc.dist")

  strm <- ShortRead::FastqStreamer(fastq_file)
  while (nreads <- length(res <- ShortRead::yield(strm))) {
    match.idx <- sapply(strands, function(x) {
      data.frame(matrix(
        ncol = length(stat.names), 
        nrow = nreads, 
        dimnames = list(NULL,stat.names)
      ))
    }, simplify = FALSE)

    for (strand in sort(strands, decreasing = TRUE)) {
      if (strand == "-") {
        res <- Biostrings::reverseComplement(res)
      }

      # save valid match results to m
      pattern <- paste0(c("(", paste0(structure, collapse = ")("), ")"), collapse = "")
      m <- aregexec(pattern, ShortRead::sread(res), max.distance = pattern.max.distance)
      match.idx[[strand]]$matched <- sapply(m, function(x) {
        x[1] != -1 # valid match
      })
      res.matched <- res[match.idx[[strand]]$matched]
      m <- m[match.idx[[strand]]$matched]
      if (!length(m)) {
        next
      }

      if ("umi" %in% names(structure)) {
        umi <- sapply(regmatches(ShortRead::sread(res.matched), m), function(x) {
          x[which("umi" == names(structure)) + 1]
        })
        stopifnot(length(m) == length(umi)) # delete
      } else {
        umi <- rep("", length(m))
      }

      if ("barcode" %in% names(structure)) {
        bc <- sapply(regmatches(ShortRead::sread(res.matched), m), function(x) {
          x[which("barcode" == names(structure)) + 1]
        })
        bc.exact <- bc %in% bc.list
        if (bc.max.distance >= 1 && length(bc[!bc.exact])) {
          bc.dists <- sapply(bc[!bc.exact], function(x) {
            adist(x, bc.list)
          })
          bc.fuzzy <- apply(bc.dists, 2, function(x) {
            ifelse(min(x) < bc.max.distance && sum(x == min(x)) == 1, which.min(x),
              -1
            )
          })
          bc[!bc.exact][bc.fuzzy != -1] <- bc.list[bc.fuzzy[bc.fuzzy != -1]]

          bc.matched <- bc.exact
          bc.matched[!bc.exact][bc.fuzzy != -1] <- TRUE
        } else {
          bc.matched <- bc.exact
        }
         match.idx[[strand]]$matched[match.idx[[strand]]$matched][!bc.matched] <- FALSE
        readseq <- Biostrings::DNAStringSet(sapply(
          regmatches(ShortRead::sread(res.matched), m),
          function(x) {
            x[which("seq" == names(structure)) + 1]
          }
        ))[bc.matched]
        readid <- Biostrings::BStringSet(paste0(bc, "_", umi, "#", ShortRead::id(res.matched)))[bc.matched]
        readqual <- sapply(regmatches(Biostrings::quality(res.matched), m), function(x) {
          x[which("seq" == names(structure)) + 1]
        })[bc.matched]
        readqual <- ShortRead::SFastqQuality(readqual)
      } else {
        readseq <- Biostrings::DNAStringSet(sapply(
          regmatches(ShortRead::sread(res.matched), m),
          function(x) {
            x[which("seq" == names(structure)) + 1]
          }
        ))
        readid <- Biostrings::BStringSet(paste0(umi, "#", ShortRead::id(res.matched)))
        readqual <- sapply(regmatches(Biostrings::quality(res.matched), m), function(x) {
          x[which("seq" == names(structure)) + 1]
        })
        readqual <- ShortRead::SFastqQuality(readqual)
      }

      if (!missing(trim.after) && !is.null(trim.after)) {
        alig <- Biostrings::pairwiseAlignment(subject = trim.after, pattern = readseq,
          type = "overlap")
        readseq <- XVector::subseq(readseq, start = ifelse(BiocGenerics::score(alig) >
          10, end(pattern(alig)), 1))
        readqual <- ShortRead::SFastqQuality(XVector::subseq(readqual, start = ifelse(BiocGenerics::score(alig) >
          10, end(pattern(alig)), 1)))
      }

      demultiplexed_res <- new("ShortReadQ", sread = readseq, quality = readqual,
        id = readid)
      ShortRead::writeFastq(object = demultiplexed_res[BiocGenerics::width(demultiplexed_res) >
        10], mode = "a", file = out.file)
    }
#    if (length(strands) > 1) {
#      cat(as.character(ShortRead::id(res[match.idx[["+"]] & match.idx[["-"]]])),
#        sep = "\n")
#      cat("\n")
#    }
  }
  close(strm)

#  if (length(strands) > 1) {
#    stats["n.mixed.strands"] <- sum(match.idx[["+"]] & match.idx[["-"]])
#  }
#  for (strand in strands) {
#    stats[paste0("n.demultiplexed.", strand)] <- sum(match.idx[[strand]])
#  }
#  return(stats)
}

# structure = list( adaptor.p5 = 'CTACAC', barcode = '.{16}',
# adaptor.read1nspacer = 'CGCGTC.{28}GAGACAG', seq = '.*' ),

agrep_demultiplex <- function(fastq_file, bc.list, out.file = "out.fq.gz", structure = list(partial.truseqread1 = "CTACACGACGCTCTTCCGATCT",
  barcode = ".{16}", umi = ".{10}", polyT = "T{3,20}", seq = ".*"), trim.after = NULL,
  strands = c("-", "+"), pattern.max.distance = 3, bc.max.distance = 2) {
  if (length(fastq_file) == 1) {
    return(agrep_demultiplex_(fastq_file, bc.list, out.file, structure, trim.after,
      strands, pattern.max.distance, bc.max.distance))
  } else {
    bcs <- data.frame()
    for (i in fastq_file) {
      bcs <- rbind(bcs, agrep_demultiplex_(i, bc.list, out.file, structure,
        trim.after, strands, pattern.max.distance, bc.max.distance))
    }
    return(bcs)
  }
}

if (FALSE) {
  barcodes.atac.tsv <- read.csv("./barcodes.atac.tsv", header = F)$V1
  cat(format(Sys.time(), "%a %b %d %X"), "\n")
  x <- agrep_demultiplex(list.files(c("../fastq_pass", "../fastq_fail"), full.names = TRUE,
    pattern = "*.fastq.gz"), barcodes.atac.tsv, bc.max.distance = 1, pattern.max.distance = 3)
  cat(format(Sys.time(), "%a %b %d %X"), "\n")

  saveRDS(x, "stats.rds")
}
