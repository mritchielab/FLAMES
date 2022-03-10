import subprocess
import os
import sys
#import pysam

# depreciated conversion using paftools.
# new version (gtf_to_bed in gff3_to_bed module) uses custom script


def gff3_to_bed12(mm2_prog_path, gff3_file, bed12_file):
    if os.path.isfile(os.path.join(mm2_prog_path, "k8")):
        cmd = "{_k8} {_paftools} gff2bed {_gff3} > {_bed}".format(
            _k8=os.path.join(mm2_prog_path, "k8"),
            _paftools=os.path.join(mm2_prog_path, "paftools.js"),
            _gff3=gff3_file,
            _bed=bed12_file)
    elif os.path.isfile(os.path.join(os.path.dirname(mm2_prog_path), "k8")):
        cmd = "{_k8} {_paftools} gff2bed {_gff3} > {_bed}".format(
            _k8=os.path.join(os.path.dirname(mm2_prog_path), "k8"),
            _paftools=os.path.join(os.path.dirname(
                mm2_prog_path), "paftools.js"),
            _gff3=gff3_file,
            _bed=bed12_file)
    else:
        print("k8 not found in minimap2_dir or its parent folder")
        cmd = "paftools.js gff2bed {_gff3} > {_bed}".format(
            _gff3=gff3_file,
            _bed=bed12_file)
    print(cmd)
    print(subprocess.check_output([cmd], shell=True, stderr=subprocess.STDOUT))
    return


def minimap2_align(mm2_prog_path, fa_file, fq_in, sam_out, no_flank=False, bed12_junc=None):
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
        no_flank = "--splice-flank=no"
    else:
        no_flank = ""
    # align_cmd = "{_prog} -ax splice -t 12 {_others} -k14 --secondary=no {_index} {_fq} | samtools view -bS -@ 4 -m 2G -o {_out} -  ".format(\
    align_cmd = "{_prog} -ax splice -t 12 {_others} -k14 --secondary=no {_index} {_fq} -o {_out}".format(
        _prog=os.path.join(mm2_prog_path, "minimap2"), _index=fa_file, _fq=fq_in, _out=sam_out, _others=" ".join([junc_cmd, no_flank]))
    print(subprocess.check_output([align_cmd],
          shell=True, stderr=subprocess.STDOUT))
    sys.stdout.flush()


def minimap2_tr_align(mm2_prog_path, fa_file, fq_in, sam_out):
    """
    minimap2 align to transcript
    """
    align_cmd = "{_prog} -ax map-ont -p 0.9 --end-bonus 10 -N 3 -t 12 {_index} {_fq} -o {_out}".format(
        _prog=os.path.join(mm2_prog_path, "minimap2"), _index=fa_file, _fq=fq_in, _out=sam_out)
    # print align_cmd
    print(subprocess.check_output([align_cmd],
          shell=True, stderr=subprocess.STDOUT))
    sys.stdout.flush()


def check_minimap2_available(mm2_prog_path):
    check_cmd = "{_prog} --help".format(
        _prog=os.path.join(mm2_prog_path, "minimap2"))
    try:
        subprocess.check_output([check_cmd], shell=True,
                                stderr=subprocess.STDOUT)
        return True
    except(subprocess.CalledProcessError):
        return False


if __name__ == "__main__":
    mm2_path = "/Users/voogd.o/Documents/GitHub/minimap2"
    print("\tminimap2 available at real location? ",
          check_minimap2_available(mm2_path))

    print("\tminimap2 available at bad location? ",
          check_minimap2_available("/Users/voogd.o/Documents/GitHub"))
