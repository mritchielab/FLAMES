from parse_gene_anno import parse_gff_tree
from sc_longread import blocks_to_junctions, remove_similar_tr, get_gene_blocks, get_gene_flat, group_bam2isoform, group_bam2isoform_multisample
from gff3_to_fa import get_transcript_seq
import sys


def find_isoform(gff3, genome_bam, isoform_gff3, tss_tes_stat, genomefa,
                 transcript_fa, downsample_ratio, config_dict, raw_splice_isoform):
    # find isoform
    print("#### Read gene annotations", flush=True)
    chr_to_gene, transcript_dict, gene_to_transcript, transcript_to_exon = parse_gff_tree(
        gff3)

    transcript_to_junctions = {tr: blocks_to_junctions(
        transcript_to_exon[tr]) for tr in transcript_to_exon}
    remove_similar_tr(gene_to_transcript, transcript_to_exon)
    gene_dict = get_gene_flat(gene_to_transcript, transcript_to_exon)
    chr_to_blocks = get_gene_blocks(gene_dict, chr_to_gene, gene_to_transcript)

    # finding isoforms are required
    print("#### find isoforms", flush=True)
    group_bam2isoform(genome_bam, isoform_gff3, tss_tes_stat, "", chr_to_blocks, gene_dict, transcript_to_junctions, transcript_dict, genomefa,
                      config=config_dict,
                      downsample_ratio=downsample_ratio, raw_gff3=raw_splice_isoform)

    # get fasta
    chr_to_gene_i, transcript_dict_i, gene_to_transcript_i, transcript_to_exon_i = parse_gff_tree(
        isoform_gff3)
    ref_dict = {"chr_to_gene": chr_to_gene, "transcript_dict": transcript_dict,
                "gene_to_transcript": gene_to_transcript, "transcript_to_exon": transcript_to_exon}
    if not config_dict["realign_parameters"]["use_annotation"]:
        ref_dict = None
    get_transcript_seq(genomefa, transcript_fa, chr_to_gene_i, transcript_dict_i,
                       gene_to_transcript_i, transcript_to_exon_i, ref_dict=ref_dict)

    sys.stdout.flush()
    return {"transcript_dict": transcript_dict, "transcript_dict_i": transcript_dict_i}


def find_isoform_multisample(gff3, genome_bams, isoform_gff3, tss_tes_stat, genomefa,
                             transcript_fa, downsample_ratio, config_dict, raw_splice_isoform):
    """
    Multisample version of `find_isoform`. 
    Requires a list of bam files (`genome_bams`) instead of only one.
    """
    # find isoform
    print("#### Read gene annotations", flush=True)
    chr_to_gene, transcript_dict, gene_to_transcript, transcript_to_exon = parse_gff_tree(
        gff3)

    transcript_to_junctions = {tr: blocks_to_junctions(
        transcript_to_exon[tr]) for tr in transcript_to_exon}
    remove_similar_tr(gene_to_transcript, transcript_to_exon)
    gene_dict = get_gene_flat(gene_to_transcript, transcript_to_exon)
    chr_to_blocks = get_gene_blocks(gene_dict, chr_to_gene, gene_to_transcript)

    # finding isoforms are required
    print("#### find isoforms", flush=True)
    group_bam2isoform_multisample(genome_bams, isoform_gff3, tss_tes_stat, "", chr_to_blocks, gene_dict, transcript_to_junctions, transcript_dict, genomefa,
                                  config=config_dict,
                                  downsample_ratio=downsample_ratio, raw_gff3=raw_splice_isoform)

    chr_to_gene_i, transcript_dict_i, gene_to_transcript_i, transcript_to_exon_i = parse_gff_tree(
        isoform_gff3)
    ref_dict = {"chr_to_gene": chr_to_gene, "transcript_dict": transcript_dict,
                "gene_to_transcript": gene_to_transcript, "transcript_to_exon": transcript_to_exon}
    if not config_dict["realign_parameters"]["use_annotation"]:
        ref_dict = None
    get_transcript_seq(genomefa, transcript_fa, chr_to_gene_i, transcript_dict_i,
                       gene_to_transcript_i, transcript_to_exon_i, ref_dict=ref_dict)

    sys.stdout.flush()
    return {"transcript_dict": transcript_dict, "transcript_dict_i": transcript_dict_i}


if __name__=="__main__":
    data  = "/Users/voogd.o/Documents/FLAMESintermediate/SIRV/"
    genome_bam = data + "FLAMESout/align2genome.bam"
    gff3 = data + "SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf"
    genomefa = data +"SIRV_isoforms_multi-fasta_170612a.fasta"
    config_file = data +  "SIRV_config_default.json"

    out = "/Users/voogd.o/Documents/FLAMESintermediate/SIRV/flamesc++out/py/"
    isoform_gff3 = out + "isoform_annotated.gff3"
    tss_tes_stat = out +"tss_tes.bedgraph"
    transcript_fa = out + "transcript_assembly.fa"
	
    config = {'comment': 'this is the default config for SIRV spike-in data. use splice annotation on alignment.', 'pipeline_parameters': {'seed': 2022, 'do_genome_alignment': True, 'do_isoform_identification': True, 'bambu_isoform_identification': False, 'do_read_realignment': True, 'do_transcript_quantification': True}, 'barcode_parameters': {'max_edit_distance': 2, 'has_UMI': False}, 'isoform_parameters': {'generate_raw_isoform': True, 'max_dist': 10, 'max_ts_dist': 100, 'max_splice_match_dist': 10, 'min_fl_exon_len': 40, 'max_site_per_splice': 3, 'min_sup_cnt': 10, 'min_cnt_pct': 0.01, 'min_sup_pct': 0.2, 'bambu_trust_reference': True, 'strand_specific': 1, 'remove_incomp_reads': 5, 'downsample_ratio': 1}, 'alignment_parameters': {'use_junctions': True, 'no_flank': True}, 'realign_parameters': {'use_annotation': True}, 'transcript_counting': {'min_tr_coverage': 0.75, 'min_read_coverage': 0.75}}

    find_isoform(gff3, genome_bam, isoform_gff3, tss_tes_stat, genomefa, transcript_fa, 1, config, "")