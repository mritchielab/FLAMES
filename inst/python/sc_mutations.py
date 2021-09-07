#!/usr/bin/env python

import os
from bam_mutations import parse_gff_tree, get_gene_flat, get_gene_blocks, get_all_SNV_table


def sc_mutations(fa_f, bam_short, out_dir, barcode_tsv, gff_f=None, known_positions=None, min_cov=100, report_pct=(0.15,0.85), test_mode=False):
    """
    fa_f: reference genome fa file (?)
    bam_short: bam file for short reads (?) 
    out_dir: FLAMES output folder containing align2genome.bam and isoform_annotated.gff3
             outputs would be saved to this folder.
    barcode_tsv: barcodes.tsv file
    known_positions: (?)

    outputs (saved to out_dir):
        ref_cnt.csv.gz
        alt_cnt.csv.gz
        allele_stat.csv.gz
        freq_summary.csv
    """

    min_cov = int(min_cov) # reticulate will convert min_cov to a float.

    # Check inputs

    # would a list of tuples be better?
    if known_positions is None:
        known_positions = ()
    else:
        assert len(known_positions) % 2 == 0, \
            "known_positions should have even length, e.g. ['chr1', 123, 'chr1', 124, 'chrX', 567]"
        assert all(isinstance(known_positions[i], (float, int)) for i in range(1, len(known_positions), 2)), \
            "Positions should be numeric value"
        known_positions = [(known_positions[i], int(known_positions[i + 1])) for i in range(0, len(known_positions), 2)]
    
    assert os.path.isfile(fa_f), "Reference genome (fa_f) not found!"
    assert (bam_short is None) or (os.path.isfile(bam_short)), "The specified short read bam file could not be found"
    assert os.path.isfile(barcode_tsv), "barcode file (barcode_tsv) not found!"
    assert len(report_pct) == 2, "report_pct should contain 2 numbers, e.g. (0.15, 0.85)"
    assert 0<=report_pct[0]<=1 and 0<=report_pct[1]<=1, "values in report_pct should be between 0 and 1, e.g. (0.15, 0.85)" 

    cb_seq_dict = dict((it.strip().split("-")[0], it.strip().split("-")[0]) for it in open(barcode_tsv))
    bam_in = os.path.join(out_dir, "align2genome.bam")

    assert os.path.isfile(bam_in), "align2genome.bam not found under " + out_dir

    if not gff_f:
        print "Using isoform_annotated.gff3 ..."
        gff_f = os.path.join(out_dir, "isoform_annotated.gff3")
        assert os.path.isfile(gff_f), "Reference annotation not provided, isoform_annotated.gff3 not found under " + out_dir
    chr_to_gene, transcript_dict, gene_to_transcript, transcript_to_exon = parse_gff_tree(gff_f)

    if test_mode:
        chr_to_gene = {chr_to_gene.keys()[0]: chr_to_gene[chr_to_gene.keys()[0]]}


    gene_dict = get_gene_flat(gene_to_transcript, transcript_to_exon)
    chr_to_blocks = get_gene_blocks(gene_dict, chr_to_gene, gene_to_transcript)
    get_all_SNV_table(bam_in, chr_to_blocks, transcript_to_exon, fa_f, out_dir, cb_seq_dict, bam_short, known_positions, min_cov, report_pct)

    return None


if __name__ == '__main__':
    """
    For debugging only
    """
    sc_mutations("/Volumes/Mattlab/LuyiTian/Index/GRCh38.primary_assembly.genome.fa", None,
                 "/Volumes/Mattlab/Changqing", "/Volumes/Mattlab/Changqing/barcodes.tsv", ["chr18", 63318364], True)
