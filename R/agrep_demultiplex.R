#fq_files <- list.files(list.files(path = folder, include.dirs=FALSE, pattern = "\\.(fastq|fq)(\\.gz)?$"))
agrep_demultiplex_ <- function(fastq_file, out.file = "out.fq.gz", pattern = '(.*)(CTCTGTC.{28}GACGCG)(.{16})(GTGTAG)', max.distance = 3, trim.after = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG") {
  bcs <- c()
  strm <-  FastqStreamer(fastq_file)
  while (length(res <- yield(strm))) {
    m <- aregexec(pattern, sread(res), max.distance = max.distance)

    res <- res[sapply(m, function(x){x[1] != -1})]
    m <- m[sapply(m, function(x){x[1] != -1})]

    bc <- BStringSet(sapply(regmatches(sread(res), m), function(x){x[4]}))
    readseq <- DNAStringSet(sapply(regmatches(sread(res), m), function(x){x[2]}))
    readid <- BStringSet(paste(bc, id(res), sep = '#'))
    
    alig <- pairwiseAlignment(subject = trim.after, pattern = readseq, type = "overlap")
    readseq <- subseq(readseq, start = ifelse(score(alig) > 10, end(pattern(alig)), 1))
    readqual <- sapply(regmatches(quality(res), m), function(x){x[2]})
    readqual <- SFastqQuality(subseq(readqual, start = ifelse(score(alig) > 10, end(pattern(alig)), 1)))

    demultiplexed_res <- new("ShortReadQ", sread = readseq, quality = readqual, id = readid)
    writeFastq(object = demultiplexed_res[width(demultiplexed_res) > 10], mode = "a", file = out.file)

    bcs <- append(bcs, bc)
  }

  return(bcs)
}

agrep_demultiplex <- function(fastq_file, out.file = "out.fq.gz", pattern = '(.*)(CTCTGTC.{28}GACGCG)(.{16})(GTGTAG)', max.distance = 3, trim.after = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG") {
  if (length(fastq_file) == 1) {
    return(agrep_demultiplex_(fastq_file, out.file, pattern, max.distance, trim.after))
  } else {
    bcs <- c()
    for (i in fastq_file) {
    bcs <- append(bcs, agrep_demultiplex_(i, out.file, pattern, max.distance, trim.after))
    }
    return(bcs)
  }
}

x <- agrep_demultiplex(
  sample(list.files("../pass/fastq_pass", full.names = TRUE)[grepl("\\.fastq\\.gz$", list.files("../pass/fastq_pass"))], size = 100)
)
