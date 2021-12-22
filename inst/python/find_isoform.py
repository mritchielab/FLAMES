from parse_gene_anno import parse_gff_tree
from sc_longread import blocks_to_junctions, remove_similar_tr, get_gene_blocks, get_gene_flat, group_bam2isoform
from gff3_to_fa import get_transcript_seq
from parse_config import parse_json_config

def log_params(chr_to_blocks, gene_dict, transcript_to_junctions, transcript_dict, outfile="group_bam2isoform_params_py.txt"):
    with open(outfile, "w") as file:
        file.write("chr_to_blocks:\n")
        for chrom in chr_to_blocks:
            file.write("\t{}:\n".format(chrom))
            for block in chr_to_blocks[chrom]:
                file.write("\t\t({},{})\n".format(block.s, block.e))

        file.write("\ngene_dict:\n")
        for gene in gene_dict:
            file.write("\t{}:\n".format(gene, gene_dict[gene]))
            for entry in gene_dict[gene]:
                file.write("\t\t({},{})\n".format(entry[0],entry[1]))

        file.write("\ntranscript_to_junctions:\n")
        for transcript in transcript_to_junctions:
            file.write("\t{}:{}\n".format(transcript, transcript_to_junctions[transcript]))

        file.write("\ntranscript_dict:\n")
        for tr in transcript_dict:
            file.write("\t{}:\n".format(tr))
            file.write("\t\t({},{})\n".format(transcript_dict[tr].start, transcript_dict[tr].end))

def log_gff_data(chr_to_gene, transcript_dict, gene_to_transcript, transcript_to_exon, filename):
    """
        logs (to filename) the contents of everything that has been extracted from a GFF file
        purely for debugging purposes
    """
    
    with open(filename, "w") as file:
        file.write("chr_to_gene: (size {})\n".format(len(chr_to_gene)))
        for ch in chr_to_gene:
            file.write("\t{}: (size {}) {}\n".format(ch, len(chr_to_gene[ch]), chr_to_gene[ch]))

        file.write("\ntranscript_dict: (size {})\n".format(len(transcript_dict)))
        for tr in transcript_dict:
            file.write("\t{}:{}\n".format(tr, transcript_dict[tr]))
        
        file.write("\ngene_to_transcript: (size {})\n".format(len(gene_to_transcript)))
        for gene in gene_to_transcript:
            file.write("\t{}: (size {}) {}\n".format(gene, len(gene_to_transcript[gene]), gene_to_transcript[gene]))
        
        file.write("\ntranscript_to_exon: (size {})\n".format(len(transcript_to_exon)))
        for tr in transcript_to_exon:
            file.write("\t{}: (size {}) {}\n".format(tr, len(transcript_to_exon[tr]), transcript_to_exon[tr]))

def find_isoform(gff3, genome_bam, isoform_gff3, tss_tes_stat, genomefa, 
                transcript_fa, downsample_ratio, config_dict, raw_splice_isoform):
    # find isoform
    print "#### Read genne annotations"
    chr_to_gene, transcript_dict, gene_to_transcript, transcript_to_exon = parse_gff_tree(gff3)
    log_gff_data(chr_to_gene, transcript_dict, gene_to_transcript, transcript_to_exon, "initial_gff_log_py.txt")

    print "len(chr_to_gene):"
    print len(chr_to_gene)
    for chrom in chr_to_gene:
        print "\t{}:{}\n".format(chrom, len(chr_to_gene[chrom]))
    print "len(transcript_dict):"
    print len(transcript_dict)
    print "len(gene_to_transcript):"
    print len(gene_to_transcript)
    for gene in gene_to_transcript:
        print "\t{}:{}\n".format(gene, len(gene_to_transcript[gene]))
    print "len(transcript_to_exon):"
    print len(transcript_to_exon)
    
    transcript_to_junctions = {tr: blocks_to_junctions(transcript_to_exon[tr]) for tr in transcript_to_exon}
    print "len(transcript_to_junctions):"
    print len(transcript_to_junctions)
    for transcript in transcript_to_junctions:
        print "\t{}:{}\n".format(transcript, transcript_to_junctions[transcript])

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
    log_params(chr_to_blocks, gene_dict, transcript_to_junctions, transcript_dict)
    group_bam2isoform(genome_bam, isoform_gff3, tss_tes_stat, "", chr_to_blocks, gene_dict, transcript_to_junctions, transcript_dict, genomefa,
    config=config_dict["isoform_parameters"], 
    downsample_ratio=downsample_ratio, raw_gff3=None)
    #raw_gff3=raw_splice_isoform if config_dict["global_parameters"]["generate_raw_isoform"] else None)


    # get fasta
    #print "### generate transcript fasta file", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    chr_to_gene_i, transcript_dict_i, gene_to_transcript_i, transcript_to_exon_i = parse_gff_tree(isoform_gff3)
    log_gff_data(chr_to_gene_i, transcript_dict_i, gene_to_transcript_i, transcript_to_exon_i, "final_gff_log_py.txt")

    ref_dict = {"chr_to_gene":chr_to_gene, "transcript_dict":transcript_dict, "gene_to_transcript":gene_to_transcript, "transcript_to_exon":transcript_to_exon}
    if not config_dict["realign_parameters"]["use_annotation"]:
        ref_dict = None
    get_transcript_seq(genomefa, transcript_fa, chr_to_gene_i, transcript_dict_i, gene_to_transcript_i, transcript_to_exon_i,ref_dict=ref_dict)

    with open("isoform_objects_py.txt", 'w') as file:
        file.write("transcript_dict:\n")
        for entry in transcript_dict:
            val = transcript_dict[entry]
            file.write("\t{}:{},{},{}\n".format(entry, val.chr, val.start, val.end))
        file.write("transcript_dict_i:\n")
        for entry in transcript_dict_i:
            val = transcript_dict_i[entry]
            file.write("\t{}:{},{},{}\n".format(entry, val.chr, val.start, val.end))
        
    return {"transcript_dict": transcript_dict, "transcript_dict_i": transcript_dict_i}
