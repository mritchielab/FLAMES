import bamnostic as bs

def bam_read(bam_in, ch, s, e):
    print "started bam_read"

    # read a bamfile
    bamfile = bs.AlignmentFile(bam_in, "rb")

    it_region = bamfile.fetch(ch, s, e)

    recnum = 0
    for rec in it_region:
        print "looking at rec {}".format(recnum)
        recnum = recnum + 1

        print "\t\t{}: ({},{})".format(rec.reference_name, rec.reference_start, rec.reference_start + rec.reference_length)