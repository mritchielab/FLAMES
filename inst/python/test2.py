import pysam
from parse_gene_anno import parse_gff_tree
from sc_longread import blocks_to_junctions

def add(x, y):
    return x + y

def t(gff3):
    chr_to_gene, transcript_dict, gene_to_transcript, transcript_to_exon = parse_gff_tree(gff3)
    # so chr_to_gene is literally just a dictionary of labels with the values being the same label? 
    transcript_to_junctions = {tr: blocks_to_junctions(transcript_to_exon[tr]) for tr in transcript_to_exon}
    print(transcript_to_junctions)
    print("hello\n")
    print(type(transcript_to_junctions))
    return transcript_to_junctions

def t2(transcript_to_junctions):
    print(transcript_to_junctions)

def t3(first, lskwargs, **kwargs):
    print first
    print "now come the kwargs"
    for key in kwargs.keys():
        print key, kwargs[key]
    print(lskwargs)
    

