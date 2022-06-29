# flames_test_function_r <- function() {
#     library(Rcpp)
#     cat("calling test function\n")
    
# }

# @export
test_find_isoform <- function() {
    x <- test_fi()
    # TODO: finish some tests for this. See if it can be integrated into testthat???

    # now we need to test our output files against python version 
    # (isoform_gff3, transcript_fa, raw_splice, as well as find_isoform return)
    if (length(x$isoform$transcript_dict) != 69)
        cat("transcript dict is the wrong size (need 69, have ",
            length(x$isoform$transcript_dict), ")\n")
    if (length(x$isoform$transcript_dict_iso) != 43)
        cat("transcript_dict_iso is the wrong size (need 43, have ",
            length(x$isoform$transcript_dict_iso), ")\n")

    isoform_gff3 <-
        read.csv(x$isoform_gff3, sep="\t", header=FALSE, 
                col.names=c("sequence", "source", "feature", "start", "end",
                            "score", "strand", "phase", "attributes"),
                            comment.char = "#")


    # what are we actually going to test here?
    print(x$transcript_fa)
}