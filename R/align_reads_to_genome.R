#' Title
#'
#' DESC
#' 
#' @param name desc
#'
#' @param name desc
#' @export
align_reads_to_genome <- function(gene_anno, input_fq, genome_fa, genome_bam, minimap2_dir, out_dir, use_junctions, no_flank) {
    print ("### align reads to genome using minimap2") # need time?
    bed <- paste(out_dir, "tmp.splice_anno.bed12", sep=.Platform$file.sep)
    bam <- paste(out_dir, "tmp.align.bam", sep=.Platform$file.sep)

 #   callBasilisk(flames_env, function (gff3, infq, genomefa, minimap2_dir, tmp_bed, tmp_bam, use_junctions, no_flank) {
    if (use_junctions) {
        gff3_to_bed12(minimap2_dir, gene_anno, bed)
    }
    minimap2_align(minimap2_dir, genome_fa, input_fq, bam, no_flank=no_flank, bed12_junc= if (use_junctions) bed else NULL)
    samtools_sort_index(bam, genome_bam)
    file.remove(tmp_bam)
    if (use_junctions) file.remove(bed)
    #}, gff3=gene_anno, infq=input_fq, genomefa=genome_fa, minimap2_dir=minimap2_dir, 
    #    tmp_bed=bed, tmp_bam=bam, use_junctions=use_junctions, no_flank=no_flank)
}

