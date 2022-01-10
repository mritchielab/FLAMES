#include "GeneAnnoParser.h"

// int
// main()
// {
//     GeneAnnoParser * 
//     geneAnnoParser = new GeneAnnoParser("data/isoform_annotated.gff3");
//     geneAnnoParser->parse();
// }

/*****************************************************************************/

GeneAnnoParser::GeneAnnoParser(std::string filename)
{
    this->filename = filename;
    this->isGTF = true;
    this->annotationSource = "Ensembl";
}

/*
    parses through the entire file
*/
GFFData
GeneAnnoParser::parse()
{
    std::cout << "started parse()\n";
    // set up the file parser
    GFFParser * fileParser = new GFFParser(filename, "GTF");

    // work out which parsing function should be used
    void (GeneAnnoParser::*parseFunction)(GFFRecord*) =
        this->selectParseFunction();

    // start parsing through
    GFFRecord rec = fileParser->parseNextRecord();
    while (!(fileParser->isEmpty())) {
        if (!rec.broken) {
            (this->*parseFunction)(&rec);
        }
        rec = fileParser->parseNextRecord();
    }

    this->gffData.removeTranscriptDuplicates();
    return this->gffData;
}

/*
    selects the appropriate function to point to for parsing,
    based on whether the file is a GTF and what its annotation source is
*/
ParseFunction
GeneAnnoParser::selectParseFunction()
{
    ParseFunction parseFunction;
    if (this->isGTF) {
        parseFunction = &GeneAnnoParser::parseGTF;
    } else {
        if (this->annotationSource == "Ensembl") {
            parseFunction = &GeneAnnoParser::parseEnsembl;
        } else {
            parseFunction = &GeneAnnoParser::parseGENCODE;
        }
    }
    return parseFunction;
}

/*
    parses a single record of a GTF file
*/
void
GeneAnnoParser::parseGTF(GFFRecord * rec)
{
    std::cout << "running parseGTF on " << rec->feature << "\n";
    // first check that the record has a gene_id
    if (!rec->hasAttribute("gene_id")) {
            std::cout << "no gene_id on a " << rec->feature << "\n";
            // Rcpp::warning("Record did not have 'gene_id' attribute: %s", rec->format_attributes());
            return;
    }
    auto gene_id = rec->attributes["gene_id"];
    if (!vectorContains(gene_id, this->gffData.chr_to_gene[rec->seqname])) {
        this->gffData.chr_to_gene[rec->seqname].push_back(gene_id);
    }


    if (rec->feature == "gene") {
        // found a gene
        // skip the rest of the work
        return;
    }

    if (!rec->hasAttribute("transcript_id")) {
        // Rcpp::warning("Feature did not have 'transcript_id' attribute: %s", rec->format_attributes());
        return;
    }
    auto transcript_id = rec->attributes["transcript_id"];
    gffData.gene_to_transcript[gene_id].push_back(transcript_id);
    StartEndPair startEnd = rec->feature == "transcript" ? 
        (StartEndPair){rec->start-1, rec->end} : 
        (StartEndPair){-1, 1};
    std::cout << "added to transcript_dict\n";
    this->gffData.transcript_dict[transcript_id] = {
        rec->seqname,
        startEnd.start,
        startEnd.end,
        rec->strand,
        gene_id
    };

    if (rec->feature == "transcript") {
        // found a transcript
        // no more processing to do
        return;
    }

    // found an exon
    std::cout << "added to transcript_to_exon\n";
    this->gffData.transcript_to_exon[transcript_id].push_back({
        rec->start - 1,
        rec->end
    });
}

/*
    parses a single record of an Ensembl GFF file
*/
void
GeneAnnoParser::parseEnsembl(GFFRecord * rec)
{
    std::cout << "running parseEnsembl on " << rec->feature << "\n";

    if (rec->feature == "gene") {
        // found a gene
        this->gffData.chr_to_gene[rec->seqname].push_back(rec->attributes["gene_id"]);
    } else if (rec->feature == "transcript") {
        // found a transcript
    } else if (rec->feature == "exon") {
        // found an exon
    }
    auto gene_id = rec->attributes["gene_id"];
    auto transcript_id = rec->attributes["transcript_id"];

    if (rec->feature == "transcript") {
        auto parent = rec->attributes["Parent"];

        this->gffData.gene_to_transcript[gene_id].push_back(transcript_id);
    }
}

/*
    parses a single record of a GENCODE GFF file
*/
void
GeneAnnoParser::parseGENCODE(GFFRecord * rec)
{
    std::cout << "running parseGENCODE\n";
}

/*****************************************************************************/
