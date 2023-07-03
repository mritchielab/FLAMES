tmp <- function(demultiplex_stats, flagstat, sce_transcripts) {
sce_transcripts <- scuttle::addPerCellQC(sce_transcripts)
scater::plotColData(scm_lib20_transcripts, y = 'detected')
scater::plotColData(scm_lib20_transcripts, y = 'sum')
}
