def trimgff3(input, out):
    f_in = open(input, 'r')
    f_out = open(out, 'w')
    #t_i = 0
    for line in f_in:
        if (line[0]) == "#":
            f_out.write(line)
            continue
        split_line = line.split('\t')
        #print(split_line)
        if split_line[0] == "chr22" or split_line[0] == "chr1":
            print('here')
            f_out.write('\t'.join(split_line))

        #if (t_i > 10):
        #    break
        #t_i += 1
    f_in.close()
    f_out.close()

def trimfa(input, out):
    # trim the fa to chr22 and the chr after it
    f_in = open(input, 'r')
    f_out = open(out, 'w')
    chr_count = 0
    writting = False
    for line in f_in:
        if (line[0] == ">"):
            if (writting):
                writting = False
            if ("chr22" in line or "chr1" in line):
                writting = True
        if (writting):
            f_out.write(line)
    f_in.close()
    f_out.close()

        

if __name__=="__main__":
    input = "/Users/voogd.o/Documents/gencode.v33.annotation.gff3"
    output = "/Users/voogd.o/Documents/FlamesR/inst/data/genocodeshortened.v33.annotation.gff3"
    trimgff3(input, output)

    input = "/Users/voogd.o/Documents/GRCh38.primary_assembly.genome.fa"
    output = "/Users/voogd.o/Documents/FlamesR/inst/data/GRCh38short.primary_assembly.genome.fa"
    trimfa(input, output)


    
