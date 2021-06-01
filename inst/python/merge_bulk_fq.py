# create merged fastq file for bulk long read RNAseq
from collections import Counter
import gzip
import os
import random
import string

class Fastq(object):
    def __init__(self, args):
        self.name = args[0][1:]
        self.seq = args[1]
        self.qual = args[3]

    def __repr__(self):
        return "Fastq({name})".format(name=self.name)

    def __str__(self):
        return "@{name}\n{seq}\n+\n{qual}\n".format(
            name=self.name, seq=self.seq, qual=self.qual)
    def to_fa(self):
        return ">{name}\n{seq}\n".format(name=self.name, seq=self.seq)


def readfq(fq):
    '''
    Enemerator over a fastq file to generate Fastq objects containing the relevant fastq
    entries (header, sequence and quality) of a sequence.
    Fastq objects are yielded one at a time
    '''
    record = []
    openFunc = gzip.open if fq.endswith(".gz") else open
    for line in openFunc(fq):
        record.append(line.strip())
        if len(record) == 4:
            yield Fastq(record)
            record = []


def merge_bulk_fq(fq_dir, out_fq):
    """
    Merges all fastq files in the fq_dir into a single file - out_fq.
    For all fastq files in the given fq_dir, read in each sequence, record where it came from
    and store the sequence data in out_fq. anno_csv stores the corresponding pseudo barcode
    for each file, whereas out_fq stores the sequences in the form:
        @name (pseudo_bc_NNN#headerline)\n, sequence\n+\n, quality\n
    """
    # get the fastq file names
    fq_names = [it for it in os.listdir(fq_dir) if ("fq" in it) or ("fastq" in it)]
    fq_dict = {}
    # populate the fastq dictionary with the absolute paths
    for f in fq_names:
        fq_dict[f] = os.path.join(fq_dir,f)
        
    merged_fq = gzip.open(out_fq, 'wb')
    #pseudo_bc_dict = {}
    fq_cnt = Counter()
    #random.seed(2333666)

    ## double comment lines indicate ways to use less time

    ## anno_file = open(anno_csv, "w")
    ## anno_file.write("file_name,pseudo_barcode\n")
    #for i in range(len(fq_dict.keys)):
    for a_fq in fq_dict:
        # generate random barcodes to represent each fastq file
        # this section can be sped up by producing a random barcode, and then
        # simply incrementing it? as opposed to making a random one each time, 
        # and checking if it isn't already in our barcode dict
        #while True:
        #    pseudo_bc = ''.join(random.choice(string.ascii_uppercase) for _ in range(16))
        #    if pseudo_bc not in pseudo_bc_dict:
        #        print a_fq, pseudo_bc
        #        break
        # populate the dictionary with the new random fastq key
        #pseudo_bc_dict[pseudo_bc] = a_fq
        ## anno_file.write("{},{}\n".format(a_fq, pseudo_bc))
        # for each sequence in the current fastq file, write it the merged_fq file
        for rec in readfq(fq_dict[a_fq]):
            rec.name = "{}_NNN#{}".format(a_fq, rec.name)
            merged_fq.write(rec.__str__())
            fq_cnt[a_fq] += 1
            
    merged_fq.close()
    ## anno_file.close()
    # write the file names and pseudo barcodes again into the anno file
    ## this whole loop can be deleted if ## comments are shown
    #with open(anno_csv,"w") as f:
    #    f.write("file_name,pseudo_barcode\n")
    #    for pseudo_bc in pseudo_bc_dict:
    #        f.write("{},{}\n".format(pseudo_bc_dict[pseudo_bc],pseudo_bc))

    # print the counts for each fastq file, indicating the number of sequences in each
    for i in fq_cnt:
        print i, fq_cnt[i]
