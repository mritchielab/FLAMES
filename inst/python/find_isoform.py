from parse_gene_anno import parse_gff_tree
from sc_longread import blocks_to_junctions, remove_similar_tr, get_gene_blocks, get_gene_flat, group_bam2isoform
from gff3_to_fa import get_transcript_seq
from parse_config import parse_json_config

def find_isoform(gff3, genome_bam, isoform_gff3, tss_tes_stat, genomefa, 
                transcript_fa, downsample_ratio, config_dict, raw_splice_isoform):
    # find isoform
    print "#### Read genne annotations"
    chr_to_gene, transcript_dict, gene_to_transcript, transcript_to_exon = parse_gff_tree(gff3)

    print "len(chr_to_gene):"
    print len(chr_to_gene)
    for chrom in chr_to_gene:
        print "\t{}:{}".format(chrom, len(chr_to_gene[chrom]))
    print "len(transcript_dict):"
    print len(transcript_dict)
    print "len(gene_to_transcript):"
    print len(gene_to_transcript)
    for gene in gene_to_transcript:
        print "\t{}:{}".format(gene, len(gene_to_transcript[gene]))
    print "len(transcript_to_exon):"
    print len(transcript_to_exon)
    
    transcript_to_junctions = {tr: blocks_to_junctions(transcript_to_exon[tr]) for tr in transcript_to_exon}
    print "len(transcript_to_junctions):"
    print len(transcript_to_junctions)
    remove_similar_tr(gene_to_transcript, transcript_to_exon)
    gene_dict = get_gene_flat(gene_to_transcript, transcript_to_exon)

    print "len(gene_dict):"
    print len(gene_dict)
    chr_to_blocks = get_gene_blocks(gene_dict, chr_to_gene, gene_to_transcript)
    print "len(chr_to_blocks):"
    print len(chr_to_blocks)
    for chrom in chr_to_blocks:
        print "\tlen(chr_to_blocks[{}]: {}".format(chrom, len(chr_to_blocks[chrom]))
    
    # finding isoforms are required
    print "#### find isoforms"
    group_bam2isoform(genome_bam, isoform_gff3, tss_tes_stat, "", chr_to_blocks, gene_dict, transcript_to_junctions, transcript_dict, genomefa,
    config=config_dict["isoform_parameters"], 
    downsample_ratio=downsample_ratio, raw_gff3=None)
    #raw_gff3=raw_splice_isoform if config_dict["global_parameters"]["generate_raw_isoform"] else None)


    # get fasta
    #print "### generate transcript fasta file", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    chr_to_gene_i, transcript_dict_i, gene_to_transcript_i, transcript_to_exon_i = parse_gff_tree(isoform_gff3)
    ref_dict = {"chr_to_gene":chr_to_gene, "transcript_dict":transcript_dict, "gene_to_transcript":gene_to_transcript, "transcript_to_exon":transcript_to_exon}
    if not config_dict["realign_parameters"]["use_annotation"]:
        ref_dict = None
    get_transcript_seq(genomefa, transcript_fa, chr_to_gene_i, transcript_dict_i, gene_to_transcript_i, transcript_to_exon_i,ref_dict=ref_dict)

    # with open("isoform_objects_py.txt", 'w') as file:
    #     file.write("transcript_dict:\n")
    #     file.write(transcript_dict)
    #     file.write("transcript_dict_i:\n")
    #     file.write(transcript_dict_i)

    return {"transcript_dict": transcript_dict, "transcript_dict_i": transcript_dict_i}
