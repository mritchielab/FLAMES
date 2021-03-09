
# FlamesR

## Documentation
Currently, only documentation for ```bulk_long_pipeline()``` and related functions
are needed. ```sc_long_pipeline()``` is not yet finished.
R function wrappers needing documentation are contained in each top level
list item of the below list. The python functions of the same name called by 
these wrappers are available in files of the same name in inst/python. 
Each list item contains the functions needing documentation as well as a list 
of the specific items needing documentation.

To see the functions used in context of the bulk long pipeline, the file 
long_pipeline.R contains the function ```generic_long_pipeline()``` which is
run by ```bulk_long_pipeline()``` and ```sc_long_pipeline()```.

### Documentation Needed
1. sc_long_pipeline.R
    1. sc_long_pipeline (these arguments are inherited for bulk_long_pipeline)
        1. Arguments
        1. @description
1. sc_longread_functions.R
    1. group_bam2isoform
        1. Arguments
        1. @description
        1. config @details
    1. get_gene_blocks
        1. @description
        1. Arguments
        1. @return
    1. get_gene_flat
        1. @description
        1. Arguments
        1. @return 
    1. remove_similar_tr
        1. Arguments
        1. @details
    1. blocks_to_junctions
        1. @description
        1. @details
1. parse_gene_anno
    1. parse_gff_tree
        1. @return
1. gff3_to_fa
    1. get_transcript_seq
        1.Arguments
        1. @description
1. count_tr
    1. parse-realigned_bam
        1. Arguments
        1. @description
    1. wrt_tr_to_csv
        1, Arguments
        1. @description
        1. @return
1. filter_gff
    1. annotate_filter_gff
        1.Arguments
        1. Title
        1. @description
1. src/match_cell_barcode.cpp
    1. match_cell_barcode
        1. @description
        1. @return
        1. Arguments
