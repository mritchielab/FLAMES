import subprocess
import os
#import pysam

# DOCUMENT THIS?
# edited version of gff3 to bed12
def gff3_to_bed12(mm2_prog_path, gff3_file, bed12_file):
    if mm2_prog_path != "":
        cmd = "{_k8} {_paftools} gff2bed {_gff3} > {_bed}".format(
            _k8=os.path.join(mm2_prog_path, "k8"),
            _paftools=os.path.join(mm2_prog_path, "paftools.js"),
            _gff3=gff3_file,
            _bed=bed12_file)
    else:
        cmd = "paftools.js gff2bed {_gff3} > {_bed}".format(
            _gff3=gff3_file,
            _bed=bed12_file)
    print subprocess.check_output([cmd], shell=True, stderr=subprocess.STDOUT)


def minimap2_align(mm2_prog_path, fa_file, fq_in, bam_out, no_flank=False, bed12_junc=None):
    """
    minimap2 align to genome
     -ax splice -t 12 -k14 $ref $read_fq | samtools view -bS -@ 4 -m 2G -o $out_bam -
samtools sort -@ 12 -o $sorted_bam $out_bam
samtools index $sorted_bam
    """
    if bed12_junc is not None:
        junc_cmd = "--junc-bed {} --junc-bonus 1".format(bed12_junc)
    else:
        junc_cmd = ""
    if no_flank:
        no_flank="--splice-flank=no"
    else:
        no_flank=""
    align_cmd = "{_prog} -ax splice -t 12 {_others} -k14 --secondary=no {_index} {_fq} | samtools view -bS -@ 4 -m 2G -o {_out} -  ".format(\
        _prog=os.path.join(mm2_prog_path, "minimap2"), _index=fa_file, _fq=fq_in, _out=bam_out, _others=" ".join([junc_cmd, no_flank]))
    print subprocess.check_output([align_cmd], shell=True, stderr=subprocess.STDOUT)


def samtools_sort_index(bam_in, bam_out):
    cmd = "samtools sort -@ 12 -o {_sorted_bam} {_in}".format(
        _sorted_bam=bam_out,
        _in=bam_in)
    print subprocess.check_output([cmd], shell=True, stderr=subprocess.STDOUT)
    cmd = "samtools index {_sorted_bam}".format(_sorted_bam=bam_out)
    print subprocess.check_output([cmd], shell=True, stderr=subprocess.STDOUT)


def minimap2_tr_align(mm2_prog_path, fa_file, fq_in, bam_out):
    """
    minimap2 align to transcript
    """
    align_cmd = "{_prog} -ax map-ont -p 0.9 --end-bonus 10 -N 3 -t 12 {_index} {_fq} | samtools view -bS -@ 4 -m 2G -o {_out} -  ".format(\
        _prog=os.path.join(mm2_prog_path, "minimap2"), _index=fa_file, _fq=fq_in, _out=bam_out)
    # print align_cmd
    print subprocess.check_output([align_cmd], shell=True, stderr=subprocess.STDOUT)