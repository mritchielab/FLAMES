from merge_bulk_fq import Fastq, readfq
from collections import Counter
import gzip
import os
import random
import string
import csv

def demultiplex_sockeye(fq_dir, sockeye_tsv, out_fq):

    id_bcumi_dict = dict()
    with open(sockeye_tsv) as sockey_file:
        tsv_reader = csv.reader(sockey_file, delimiter="\t") 
        
        header = next(tsv_reader)
        read_id = header.index('read_id')
        barcode = header.index('barcode')
        umi = header.index('umi')
        for record in tsv_reader:
            id_bcumi_dict[record[read_id].split("_")[0]] = (record[barcode], record[umi])

    fq_names = [it for it in os.listdir(
        fq_dir) if ("fq" in it) or ("fastq" in it)]
    fq_dict = {}
    for f in fq_names:
        fq_dict[f] = os.path.join(fq_dir, f)
    merged_fq = gzip.open(out_fq, 'wt')
    fq_cnt = Counter()
    for a_fq in fq_dict:
        for rec in readfq(fq_dict[a_fq]):
            try:
                id = rec.name.split(" ")[0]
                bc, umi = id_bcumi_dict[id]
            except KeyError as ke:
                #print(f"Could not find {id}'s barcode and UMI")
                continue
            rec.name = "{}_{}#{}".format(bc, umi, id)
            merged_fq.write(rec.__str__())
            fq_cnt[a_fq] += 1
    merged_fq.close()
    for i in fq_cnt:
        print(i, fq_cnt[i])

if __name__ == "__main__":
    demultiplex_sockeye("./sockeye_demultiplex_tests/in_fq", "./sockeye_demultiplex_tests/cell_umi_gene.tsv", "./sockeye_demultiplex_tests/out.fq.gz")
