# quantify gene
import pysam
import sys
import os
import numpy as np
#import gtfparse
import pandas as pd
from parse_gene_anno import parseGFF3
import re
from tqdm import tqdm
from collections import Counter
import fast_edit_distance
from concurrent.futures import ThreadPoolExecutor, as_completed
import multiprocessing as mp
import matplotlib.pyplot as plt

def get_read_to_gene_assignment(in_bam, in_gtf, methods):
    """
    Get gene counts from a bam file and a gtf file.
    Input:
        in_bam: bam file path
        in_gtf: gtf file path
        methods: demultiplexing methods, 'flexiplex' or 'blaze'
    Process:
        Step 1: build index 
            chr -> gene -> pos
        Step 2: Assign read to gene: read_id -> chr, pos, strand, gene_id
    Output:
        gene_idx_df: a dataframe of gene index from step 1
        read_gene_assign_df: a dataframe of read to gene assignment from step 2

    """
    # read bam file
    bam_file = pysam.AlignmentFile(in_bam, "rb")

    ## read gtf file and build index df
    gtf = parseGFF3(in_gtf)
    chr_name, gene_id, start, end = [[] for i in range(4)]
    for entry in gtf:
        if entry.type == "gene":
            chr_name.append(entry.seqid)
            gene_id.append(entry.attributes["gene_id"])
            start.append(entry.start)
            end.append(entry.end)
        else: 
            continue
    gene_idx_df = pd.DataFrame({"chr_name": chr_name, "gene_id": gene_id, 
                                "start": start, "end": end})

    # Assign read to gene based on mapping position  
    chr_names, gene_ids, bcs, umis, read_ids, positions, map_ranges = [[] for i in range(7)]
    for entry in tqdm(gene_idx_df.itertuples()):
        reads_fetch = bam_file.fetch(entry.chr_name, entry.start, entry.end)
        #
        # for each reads, get the reference_name, mapping position, strand
        for read in reads_fetch:
            if read.is_supplementary or read.is_secondary or read.is_unmapped:
                continue
            
            # get mapped position
            bc, umi, read_id, strand = flames_read_id_parser(read.query_name,methods)
            #
            # summarise the mapping position as the mapped polyT/A side of the reads
            if read.is_reverse ^ (strand == '+'):
                pos = read.reference_end
            else:
                pos = read.reference_start
            if pos > entry.end or pos < entry.start:
                continue
            else:
                map_ranges.append(read.reference_end-read.reference_start)
                bcs.append(bc)
                umis.append(umi)
                read_ids.append(read_id)
                positions.append(pos)
                chr_names.append(entry.chr_name)
                gene_ids.append(entry.gene_id)
    read_gene_assign_df = pd.DataFrame({"chr_name": chr_names, "gene_id": gene_ids, 
                            "bc": bcs, "umi": umis, "read_id": read_ids, 
                            "position": positions, "map_range": map_ranges})

    return gene_idx_df, read_gene_assign_df

def flames_read_id_parser(read_id, methods = 'flexiplex'):
    """parse the read id from FLAMES output fastq/bam file.

    Args:
        read_id (str): read id
        methods (str, optional): 'flexiplex' or 'blaze'. Defaults to 'flexiplex'.
    """
    if methods == 'flexiplex':
        # format: GGATGTTAGGTTACCT-1_AAATCAGTTCTT#de97a0c6-ff84-4528-ab10-721bc5528b57_+1of1
        bc, umi, read_id, _, _ = re.split("_|#|1of1", read_id)
        # flexiplex output is always in ployT strand of cDNA
        strand = "+"
        return bc, umi, read_id, strand
    #
    if methods == 'blaze':
        bc, umi, read_id, strand = re.split("_|#", read_id)
        return bc, umi, read_id, strand
    else:
        sys.exit("Please specify the correct methods: 'flexiplex' or 'blaze'")


def quantify_gene(in_bam, in_gtf, demulti_methods, 
                   saturation_curve_fn='', estimate_saturation=False):
    """
    Get gene counts from a bam file and a gtf file.
    Input:
        in_bam: bam file path
        in_gtf: gtf file path
        demulti_methods: demultiplexing methods, 'flexiplex' or 'blaze'
        estimate_saturation: whether to estimate saturation curve
        plot_saturation_curve: plot filename for saturation curve, if specified, 
                              `estimate_saturation` is autimatically set to True
    Output:
        gene_count_mat: a matrix of gene counts
        read_gene_assign_df: a dataframe of read to gene assignment with umi corrected
    """
    # get read to gene assignment
    gene_idx_df, read_gene_assign_df = \
        get_read_to_gene_assignment(in_bam, in_gtf, methods=demulti_methods) 
    
    # cluster reads with similar genome location (polyT side mapping position)
    read_gene_assign_df['cluster'] = \
        read_gene_assign_df.groupby(['bc', 'gene_id'])['position']\
                           .transform(map_pos_grouping)
    
    # correct umi
    read_gene_assign_df['umi_corrected'] = \
        read_gene_assign_df.groupby(['bc', 'gene_id','position'])['umi']\
                           .transform(umi_correction)

    # get gene count
    gene_count_df = \
        read_gene_assign_df.groupby(['bc', 'gene_id'])['umi_corrected'].count()
    
    # convert to matrix
    gene_count_mat = gene_count_df.reset_index().pivot(index='gene_id', 
                                                       columns='bc', 
                                                       values='umi_corrected')
    gene_count_mat = gene_count_mat.rename_axis(None, axis=0).\
                                    rename_axis(None, axis=1)

    if estimate_saturation or saturation_curve_fn:
        saturation_estimation(read_gene_assign_df.umi_corrected, 
                              saturation_curve_fn)

    return gene_count_mat, read_gene_assign_df

def map_pos_grouping(mappos, min_dist=20):
    """
    Group mapping positions into clusters. 
    Output the cluster id in the same order of input mappos.
    """
    # sort the mapping position
    mappos = np.array(mappos)

    sort_indices = np.argsort(mappos)
    mappos_sorted = mappos[sort_indices]

    # calculate the distance between adjacent mapping position
    dist = np.insert(np.diff(mappos_sorted), 0, 0)
    cluster_id = np.cumsum(dist > min_dist)
    return cluster_id[sort_indices]


def umi_correction(umis, max_ed=1):
    """
    Correct umis.
    """
    read_cnt = len(umis)
    dup_cnt = Counter(umis)
    dup_cnt = sorted(dup_cnt.most_common(), key=lambda x: (x[1], x[0]), reverse=True)
    if len(dup_cnt) == 1:
        return umis
    
    umi_mapping = {} # {putative_umi: real_umi}
    #umi_count = {} # {real_umi: count}
    for ith in range(len(dup_cnt)-1):
        umi_i, dup_i = dup_cnt[ith]
        umi_mapping[umi_i] = umi_i
        #umi_count[umi_i] = dup_i
        for jth in range(len(dup_cnt)-1, ith, -1):  # first assess the low abundant UMI
            umi_j, dup_j = dup_cnt[jth]
            if umi_j not in umi_mapping:
                if fast_edit_distance.edit_distance(umi_i, umi_j, max_ed) <= max_ed:
                    umi_mapping[umi_j] = umi_i
                    #umi_count[umi_i] += dup_j
    umi_last, dup_last = dup_cnt[-1]
    if umi_last not in umi_mapping:
        umi_mapping[umi_last] = umi_last
        #umi_count[umi_last] = dup_last

    umi_corrected = [umi_mapping[umi] for umi in umis]
    
    # return corrected umi list
    return umi_corrected

def saturation_estimation(corrected_umis, plot_fn=False, num_of_points=1000):
    """
    Estimate saturation curve.
    """
    corrected_umis = list(corrected_umis)
    # The saturation is calculated as 1-unique_umi/total_reads
    est = 1-len(set(corrected_umis))/len(corrected_umis)
    
    if not plot_fn:
        return est

    # Calculate the number of unique UMIs at each sequencing depth
    if plot_fn:
        umis = np.random.permutation(corrected_umis)
        sub_samples = np.array_split(umis, min(num_of_points, len(umis)))
        
        read_depth = np.cumsum([len(x) for x in sub_samples])
        # get the number of unique umis
        umi_counter = Counter()
        umi_counts = []
        for sub_sample in sub_samples:
            umi_counter.update(sub_sample)
            umi_counts.append(len(umi_counter))

        plt.figure(figsize=(8, 6))
        plt.plot(read_depth, umi_counts, 'o-')
        plt.xlabel("Sequencing Depth")
        plt.ylabel("Number of Unique UMIs")
        plt.title(f"Saturation Curve (Saturation = {est*100:.2f}%)")
        plt.grid(True)
        plt.savefig(plot_fn)

    return est

def _pd_parellel_transform(groupby_obj, func, num_workers=mp.cpu_count()-1, **kwargs):
    """
    Apply function to each group in a pandas groupby object in parallel.

    groupby_obj: A pandas groupby object.
    func: The function to apply to each group.
        Note: The function must take a DataFrame as its first argument and return a DataFrame.
        
    For the best performance:
        There is no garuntee that the function will be output in the same order as the groups.
    """
    # Split the DataFrame into chunks based on the group column
    chunks = [group for _, group in groupby_obj]

    # Use ThreadPoolExecutor to apply transform_chunk to each chunk in parallel
    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        futures = [executor.submit(func, chunk, **kwargs) for chunk in chunks]

        # Gather the results as they complete
        results = [future.result() for future in as_completed(futures)]

    # Concatenate the processed chunks back into a single DataFrame/Series
    result_df = pd.concat(results)


# this is the main function
def quantification(annotation, outdir, demultiplex_methods, pipeline):

    if pipeline == "sc_single_sample":
        in_bam = os.path.join(outdir, "align2genome.bam")
        out_csv = os.path.join(outdir, "gene_count.csv.gz")
        out_fig = os.path.join(outdir, "saturation_curve.png")

        gene_count_mat, _ = quantify_gene(
            in_bam, annotation, demultiplex_methods, out_fig)
        gene_count_mat.to_csv(out_csv, compression='gzip')
        return

    elif pipeline == "bulk":
        """Gene quantification is not implemented in bulk pipeline
        """
        return

    elif pipeline == "sc_multi_sample":
        try:
            assert isinstance(inbam, list)
        except AssertionError:
            raise ValueError("inbam must be a list of bam files")

        in_bams = \
            [os.path.join(outdir, f) for f in os.listdir(outdir) if f[-17:] == "_align2genome.bam"]
        for sample_bam in realign_bam:
            sys.stderr.write("parsing " + sample_bam + "...\n")
            sample = os.path.basename(sample_bam).replace('_align2genome.bam','')
            out_csv = os.path.join(outdir, sample+ "_"+"gene_count.csv.gz")
            out_fig = os.path.join(outdir, sample+ "_"+"saturation_curve.png")
            gene_count_mat, _ = quantify_gene(in_bam, annotation, out_fig)
            gene_count_mat.to_csv(out_csv, compression='gzip')
        return

    else:
        raise ValueError(f"Unknown pipeline type {pipeline}")
