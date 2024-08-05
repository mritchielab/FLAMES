import pysam
import sys
import os
import gzip
import numpy as np
import pandas as pd
from parse_gene_anno import parseGFF3
import re
from collections import Counter
import fast_edit_distance
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
import multiprocessing as mp
import matplotlib.pyplot as plt
import bisect
from collections import namedtuple

import helper
#import cProfile

def parse_gtf_to_df(in_gtf):
    """
    Parse gtf file to a dataframe.
    """
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
                                "start": start, "end": end} )
    return gene_idx_df

def get_read_to_gene_assignment(in_bam, gene_idx_df, methods):
    """
    Get gene counts from a bam file and a gtf file.
    Input:
        in_bam: bam file path
        gene_idx_df: gtf dataframe returned by parse_gtf_to_df
        methods: demultiplexing methods, 'flexiplex' or 'blaze'
    Process:
        Step 1: build index 
            chr -> gene -> pos
        Step 2: Assign read to gene: read_id -> chr, pos, strand, gene_id
    Output:
        gene_idx_df: a dataframe of gene index from step 1
        read_gene_assign_df: a dataframe of read to gene assignment from step 2

    Note: If a read is assigned to multiple genes, keep the one with the largest overlap. 
    However, if there is a tie in the overlap, discard the read. 
    """
    # read bam file
    bam_file = pysam.AlignmentFile(in_bam, "rb")

    # Assign read to gene based on mapping position  
    chr_names, gene_ids, bcs, umis, read_ids, positions_3prim, \
        positions_5prim, overlaps, read_lengths = [[] for i in range(9)]
    
    # process genes in the gene_idx_df
    for gene in gene_idx_df.itertuples():
        reads_fetch = bam_file.fetch(gene.chr_name, gene.start, gene.end)
        # for each reads, get the reference_name, mapping position, strand
        for read in reads_fetch:
            if read.is_supplementary or read.is_secondary or read.is_unmapped:
                continue
            # get mapped position
            bc, umi, read_id, strand = flames_read_id_parser(read.query_name,methods)
            
            # get the overlaps 
            if read.reference_start >= gene.start and read.reference_end <= gene.end:
                # when the read is fully mapped within the gene
                overlaps.append(read.query_alignment_length)
                if read.is_reverse ^ (strand == '+'):
                    pos3, pos5 =read.reference_end, read.reference_start
                else:
                    pos3, pos5 = read.reference_start, read.reference_end
            else:
                # when the read is not fully mapped to the gene
                ref_positions = read.get_reference_positions()# slow step
                in_gene_read_start = bisect.bisect_right(ref_positions, gene.start)
                in_gene_read_end = bisect.bisect_right(ref_positions, gene.end) -1
                
                if in_gene_read_start >= in_gene_read_end:
                    continue
                else:
                    overlaps.append(in_gene_read_end-in_gene_read_start)
                    # get the mapping position of two ends of the reads in gene
                    if read.is_reverse ^ (strand == '+'):
                        pos3, pos5 = \
                            ref_positions[in_gene_read_end], ref_positions[in_gene_read_start]
                    else:
                        pos3, pos5 = \
                            ref_positions[in_gene_read_start], ref_positions[in_gene_read_end]

            # append to the list
            positions_5prim.append(pos5)
            positions_3prim.append(pos3)
            read_lengths.append(read.reference_end-read.reference_start)
            bcs.append(bc)
            umis.append(umi)
            read_ids.append(read_id)
            chr_names.append(gene.chr_name)
            gene_ids.append(gene.gene_id)
    read_gene_assign_df = pd.DataFrame({"chr_name": chr_names, "gene_id": gene_ids,
                                        "bc": bcs, "umi": umis, "read_id": read_ids,
                                        "pos_5prim": positions_5prim, "pos_3prim": positions_3prim,
                                        "overlap": overlaps, "read_length": read_lengths})  
    # close bam file
    bam_file.close()
    # deduplication row with same gene_id and read_id
    read_gene_assign_df.drop_duplicates(subset=['gene_id', 'read_id'], inplace=True)
    # get the unambiguous read to gene assignment
    dup_mask = read_gene_assign_df.duplicated(subset='read_id', keep = False)
    unambig_df = read_gene_assign_df[~dup_mask]
    # resolve the read assigned to multipe genes
    ambig_df = read_gene_assign_df[dup_mask].sort_values(
        by=['read_id','overlap'], ascending = [True, False])
    
    # for the read assigned to multiple genes, keep the one with the largest overlap
    pre_id, pre_overlap, pre_idx = None, None, None
    row_idx_to_drop = []
    for read in ambig_df.itertuples():
        if read.read_id != pre_id:
            pre_id, pre_overlap, pre_idx = \
                read.read_id, read.overlap, read.Index
        elif read.overlap < pre_overlap:
            row_idx_to_drop.append(read.Index)
        elif read.overlap == pre_overlap:
            row_idx_to_drop.append(read.Index)
            row_idx_to_drop.append(pre_idx)
            
    recovered_ambig_df = ambig_df.drop(row_idx_to_drop)
    ## merge the unambiguous and ambiguous read to gene assignment
    read_gene_assign_df = pd.concat([unambig_df, recovered_ambig_df])


    ## sort the read_gene_assign_df
    read_gene_assign_df.sort_values(by=['chr_name', 'bc', 'gene_id', 'pos_3prim'], inplace=True)
    read_gene_assign_df[['chr_name', 'bc', 'gene_id']] =\
          read_gene_assign_df[['chr_name','bc', 'gene_id']].astype('category')

    # check if there are remaining reads assigned to multiple genes
    dup_mask = recovered_ambig_df.duplicated(subset='read_id', keep = False)
    if dup_mask.sum():
        warning_msg(f"Warning: {dup_mask.sum()} reads are assigned to multiple \
        genes. Please check the output file for details.")

    return read_gene_assign_df

def flames_read_id_parser(read_id, methods = 'flexiplex'):
    """parse the read id from FLAMES output fastq/bam file.

    Args:
        read_id (str): read id
        methods (str, optional): 'flexiplex' or 'blaze'. Defaults to 'flexiplex'.
    """
    if methods == 'flexiplex':
        # format: GGATGTTAGGTTACCT-1_AAATCAGTTCTT#de97a0c6-ff84-4528-ab10-721bc5528b57_+1of1
        bc, umi = re.split("_|#|1of1", read_id)[:2] 
        # flexiplex output is always in ployT strand of cDNA
        strand = "+"
        return bc, umi, read_id, strand
    #
    if methods == 'blaze':
        bc, umi, *_, strand = re.split("_|#", read_id)
        return bc, umi, read_id, strand
    else:
        sys.exit("Please specify the correct methods: 'flexiplex' or 'blaze'")

def quantify_gene(in_bam, in_gtf, n_process):
    """A multi-process wrapper of quantify_gene_single_process 
       Processing strategy:
        1. split the gtf file by chromosome
            2. for each chromosome, run quantify_gene_single_process
        3. combine the gene count matrix
       Input:
        in_bam: bam file path
        in_gtf: gtf file path
        n_process: number of process to run
       return:
        gene_count_mat: pd.DataFrame, a matrix of gene counts
        dup_read_lst: a list of duplicated read id
        umi_lst: a list of umi
    """
    # identify the demultiplexing methods
    bam_file = pysam.AlignmentFile(in_bam, "rb")
    first_read_id = next(bam_file).query_name
    demulti_methods = 'flexiplex' if first_read_id[-4:] == "1of1" else 'blaze'
    bam_file.close()

    # spliting the annotated gene by chrom
    in_gtf_df = parse_gtf_to_df(in_gtf)
    in_gtf_list = \
        [x for x in independent_gene_set_generator(in_gtf_df, in_bam, threads=n_process)]
    # priorities the process of the largest gene set
    in_gtf_list = sorted(in_gtf_list, key=lambda x: x.shape[0], reverse=True)

    # run quantify_gene_single_process for each chromosome
    print("Assigning reads to genes...")
    gene_count_mat_dfs, dup_read_lst, umi_lst = [], [], []
    for future in helper.multiprocessing_submit(
                            quantify_gene_single_process, 
                            in_gtf_list,
                            pbar=True,
                            pbar_unit = "gene_group",
                            preserve_order = False,
                            n_process=n_process, 
                            in_bam=in_bam, 
                            demulti_methods=demulti_methods):
        
        gene_count_mat, dup_read_lst_sub, umi_list_sub = future.result()
        gene_count_mat_dfs.append(gene_count_mat)
        dup_read_lst.extend(dup_read_lst_sub)
        umi_lst.extend(umi_list_sub)

    # combine the gene count matrix
    print("Finalising the gene count matrix ...")
    gene_count_mat = pd.concat(gene_count_mat_dfs, 
                                copy=False).fillna(0)
    
    return gene_count_mat, dup_read_lst, umi_lst

def quantify_gene_single_process(in_gtf_df, in_bam, demulti_methods, verbose=False):
    """
    Get gene counts from a bam file and a gtf file.
    Input:
        in_bam: bam file path
        in_gtf_df: gtf_df contains the subset of gene to run in a single process
        demulti_methods: demultiplexing methods, 'flexiplex' or 'blaze'
        verbose: whether to print the progress
    Output:
        gene_count_mat: a matrix of gene counts
        read_gene_assign_df: a dataframe of read to gene assignment with umi corrected
    """

    read_gene_assign_df = \
        get_read_to_gene_assignment(in_bam, in_gtf_df, methods=demulti_methods)

    # cluster reads with similar genome location (polyT side mapping position)
    if verbose:
        print("Clustering reads with similar genome locations ...")
    read_gene_assign_df.sort_values(by=['bc', 'gene_id'], inplace=True)
    cell_gene_grp = read_gene_assign_df.groupby(['bc', 'gene_id'], observed=True)

    read_gene_assign_df['cluster'] = \
        cell_gene_grp['pos_3prim'].transform(_map_pos_grouping).astype('category')
    
    # correct umi
    if verbose:
        print("Correcting UMIs ...")
    read_gene_assign_df['umi_corrected'] = \
        read_gene_assign_df.groupby(['bc','gene_id','cluster'], observed=True)['umi']\
                           .transform(_umi_correction)

    # get gene count
    if verbose:
        print("Generating per-gene UMI counts ...")
    gene_count_df = \
        read_gene_assign_df.groupby(['bc', 'gene_id'], observed=True)['umi_corrected'].nunique()
    
    # convert to matrix
    gene_count_mat = gene_count_df.reset_index().pivot(index='gene_id', 
                                                       columns='bc', 
                                                       values='umi_corrected')
    gene_count_mat = gene_count_mat.rename_axis(None, axis=0).\
                                    rename_axis(None, axis=1)

    # get list of read_id to remove
    dedup_read_lst = list_deduplicated_reads(read_gene_assign_df)

    # get list of umi (in the form of bc+umi+cluster to avoid collision)
    umi_lst = read_gene_assign_df.bc.astype(str) +\
                    read_gene_assign_df.gene_id.astype(str) + \
                    read_gene_assign_df.umi_corrected.astype(str) + \
                    read_gene_assign_df.cluster.astype(str)
    
    return gene_count_mat, dedup_read_lst, umi_lst

def _map_pos_grouping(mappos, min_dist=20):
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

def _umi_correction(umis, max_ed=1):
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

def list_deduplicated_reads(umi_corrected_df,
                            umi_col="umi_corrected", 
                            read_id_col="read_id", 
                            priority_cols="overlap", 
                            groupby_cols=["cluster", "bc", "gene_id"]):
    """
    Deduplicate reads based on umi. Keep the one with the highest overlap with the reference.

    Input:
        umi_corrected_df: A pandas DataFrame with umi_corrected column. (output of quantify_gene)
        umi_col: The column name of the umi column.
        read_col: The column name of the read column.
        priority_cols: The columns to prioritize the reads.
        groupby_cols: The columns to groupby, indicating how the umi correction has be done.
    Output:     
        list of deduplicated read id
    """
    # group by groupby_cols, and apply deduplication to each group
    # umi_deduplicated_df = umi_corrected_df.groupby(groupby_cols).apply(
    #     lambda x: x.sort_values(
    #         priority_cols, ascending=False).drop_duplicates(umi_col, keep="first")
    # )

    out_list = []
    umi_corrected_df = umi_corrected_df[groupby_cols + [umi_col, read_id_col, priority_cols]]
    for _, group in umi_corrected_df.groupby(groupby_cols + [umi_col], observed=True):
        if group.shape[0] == 1:
            continue
        else:
            # remove one read with the highest priority, randomly break ties
            read_ids = group[read_id_col].values
            priorities = group[priority_cols].values
            read_to_keep_mask= priorities==priorities.max()
            if sum(read_to_keep_mask) == 1:
                #read_to_remove = read_ids[~read_to_keep_mask]
                read_to_keep = read_ids[read_to_keep_mask]
            else:
                read_to_keep_idx = np.random.choice(np.where(read_to_keep_mask)[0])
                #read_to_remove = np.delete(read_ids, read_to_keep_idx)
                read_to_keep = read_ids[read_to_keep_idx]
            #out_list.extend(read_to_remove)
            out_list.extend(read_to_keep)

    return out_list

def saturation_estimation(corrected_umis, plot_fn=None, num_of_points=500):
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
        plt.plot(read_depth, umi_counts, '-')
        plt.xlabel("Number of reads")
        plt.ylabel("Number of Unique UMIs")
        plt.title(f"Saturation Curve (Saturation = {est*100:.2f}%)")
        plt.grid(True)
        plt.savefig(plot_fn)

    return est

def _subset_reads_from_fastq_chunk(read_chunk, read_id_set):
    """Keep reads from a chunk of fastq file
    Input:
        read_chunk: a chunk of fastq lines
        read_id_set: a set of read id to keep
    Output: 
        a list of fastq lines to keep
    """
    out_lst = []
    for i in range(0, len(read_chunk), 4):
        read_id = read_chunk[i].strip()[1:]
        if read_id in read_id_set:
            out_lst.extend(read_chunk[i:i+4])
    return out_lst

def subset_reads_from_fastq(in_fastq, out_fastq, read_id_lst, 
                            n_process,
                            chunk_size=10_000):
    """ Subset fastq file by keep only the reads in read id list
        Input:
            in_fastq: input fastq file
            out_fastq: output fastq file
            read_id_lst: list of read id to keep
            chunk_size: number of reads to read in each chunk
        Output:
            output fastq file
    """
    assert chunk_size % 4 == 0, "chunk_size must be a multiple of 4"

    if in_fastq[-3:] == ".gz":
        f_in = gzip.open(in_fastq, "rt")
    else:
        f_in = open(in_fastq, "r")
        
    if out_fastq[-3:] == ".gz":
        f_out = gzip.open(out_fastq, "wt")
    else:
        f_out = open(out_fastq, "w")

    read_chunks = helper.read_chunk_generator(f_in, chunk_size)
    read_id_lst = set(read_id_lst)
    results = helper.multiprocessing_submit( _subset_reads_from_fastq_chunk,
                                    read_chunks,
                                    n_process=12,
                                    pbar_func=lambda *x: chunk_size/4,
                                    read_id_set = read_id_lst,
                                    preserve_order=False,
                                    schduler = 'thread')
    
    for rst in results:
        f_out.writelines(rst.result())
        
    f_in.close()
    f_out.close()

    return

# this is the main function
def quantification(annotation, outdir, pipeline, 
                  n_process=12, saturation_curve=True,
                  infq=None,  sample_names=None, **kwargs):

    if pipeline == "sc_single_sample":
        if not infq:
            infq = os.path.join(outdir, "matched_reads.fastq")
        in_bam = os.path.join(outdir, "align2genome.bam")
        out_csv = os.path.join(outdir, "gene_count.csv")
        out_fig = os.path.join(outdir, "saturation_curve.png") if saturation_curve else None
        #out_read_lst = os.path.join(outdir, "duplicated_read_id.txt")
        out_fastq = os.path.join(outdir, "matched_reads_dedup.fastq")

        gene_count_mat, dedup_read_lst, umi_lst = \
                                quantify_gene(in_bam, annotation, n_process)

        #pd.DataFrame({'umi':umi_lst}).to_csv(outdir+"/umi_lst.csv")

        gene_count_mat.to_csv(out_csv)

        print("Plotting the saturation curve ...")
        saturation_estimation(umi_lst, out_fig)  

        print("Generating deduplicated fastq file ...")
        subset_reads_from_fastq(infq, out_fastq, dedup_read_lst, n_process)
        return

    elif pipeline == "bulk":
        """Gene quantification is not implemented in bulk pipeline
        """
        print("Gene quantification has not been implemented in bulk pipeline ...")
        return

    elif pipeline == "sc_multi_sample":
        if not sample_names:
            raise ValueError("Please specify the sample names")
        
        for idx, sample in enumerate(sample_names):
            sample_bam = os.path.join(outdir, sample + "_align2genome.bam")
            sys.stderr.write("parsing " + sample_bam + "...\n")
            in_fastq = infq[idx] if infq is not None else os.path.join(outdir,sample+ "_"+ "matched_reads.fastq")
            out_fastq = os.path.join(outdir,sample+ "_"+ "matched_reads_dedup.fastq")
            out_csv = os.path.join(outdir, sample+ "_"+"gene_count.csv")
            out_fig = os.path.join(outdir, sample+ "_"+"saturation_curve.png") if saturation_curve else None
            #out_read_lst = os.path.join(outdir, sample+ "_"+"deduplicated_read_id.txt")
            
            gene_count_mat, dup_read_lst, umi_lst = \
                                    quantify_gene(sample_bam, annotation, n_process)

            #pd.DataFrame({'umi':umi_lst}).to_csv(f"{outdir}/{sample}_umi_lst.csv")

            gene_count_mat.to_csv(out_csv)

            print("Plotting the saturation curve ...")
            saturation_estimation(umi_lst, out_fig)  

            print("Generating deduplicated fastq file ...")
            subset_reads_from_fastq(in_fastq, out_fastq, dup_read_lst, n_process)
        return

    else:
        raise ValueError(f"Unknown pipeline type {pipeline}")
    
def independent_gene_set_generator(gtf_df, bam_file, threads=1):
    """
    generate independent gene set from gtf_df and bam_file so that there is no
    read mapped spanning two gene sets. As a result, the gene set can be processed
    independently.
    input: 
        gtf_df: pandas.DataFrame, result from parse_gtf_to_df
        bam_file: str filename
        threads: int, number of threads to use in pysam
    """
    def _yield_genes_per_chr(gtf_df, bam):
        """
        yield genes per chromosome
        Input:
            gtf_df: pandas.DataFrame, result from parse_gtf_to_df, contains only one chromosome
            bam: pysam.AlignmentFile
        """
        genome_region = namedtuple('genome_region', ['start', 'end'])
        prev_intval = genome_region(None, None)
        yield_genes_idx = []
        for gene in gtf_df.itertuples():
            if not len(yield_genes_idx):
                prev_intval = genome_region(gene.start, gene.end)
                yield_genes_idx.append(gene.Index)
                continue
            elif gene.start < prev_intval.end:
                if next(bam.fetch(chr, gene.start,prev_intval.end), None):
                    yield_genes_idx.append(gene.Index)
                    prev_intval._replace(end = max(prev_intval.end, gene.end))
                else:
                    yield gtf_df.loc[yield_genes_idx]
                    yield_genes_idx = [gene.Index]
                    prev_intval = genome_region(gene.start, gene.end)
            elif gene.start == prev_intval.end:
                if next(bam.fetch(chr, gene.start-1,gene.start+1), None):
                    yield_genes_idx.append(gene.Index)
                    prev_intval._replace(end = gene.end)
                else:
                    yield gtf_df.loc[yield_genes_idx]
                    yield_genes_idx = [gene.Index]
                    prev_intval = genome_region(gene.start, gene.end)
            elif gene.start > prev_intval.end:
                if _check_read_spanning_region(bam, chr, prev_intval.end, gene.start+1):
                    yield_genes_idx.append(gene.Index)
                    prev_intval._replace(end = gene.end)
                else:    
                    yield gtf_df.loc[yield_genes_idx]
                    yield_genes_idx = [gene.Index]
                    prev_intval = genome_region(gene.start, gene.end)
        yield gtf_df.loc[yield_genes_idx]

    # check whether there is any read mapped fully spanning  a genome region
    def _check_read_spanning_region(bam, chr, start, end):
        """
        check whether there is any read mapped fully spanning a genome region
        (primary alignment only)
        input:
            bam: pysam.AlignmentFile
            chr: str
            start: int
            end: int
        """
        for read in bam.fetch(chr, end-1, end):
            if read.is_secondary or read.is_supplementary:
                continue
            else:
                if read.reference_start <= start:
                    return True
                else:
                    return False
    # main function
    gtf_df = gtf_df.sort_values(by=['chr_name', 'start'], ascending=True)
    bam = pysam.AlignmentFile(bam_file, 'rb', threads=threads)
    for chr in gtf_df.chr_name.unique():
        yield from _yield_genes_per_chr(gtf_df[gtf_df.chr_name==chr], bam)
    bam.close()