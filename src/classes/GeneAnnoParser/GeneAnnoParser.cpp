#include "GeneAnnoParser.h"

/*
    initialises the GeneAnnoParser
*/
GeneAnnoParser::GeneAnnoParser(std::string filename)
{
    this->filename = filename;
    // work out if the file is a GTF or not
    this->isGTF = isFileGTF(filename);
    this->annotationSource = "Ensembl";
    this->gffParser = new GFFParser(this->filename, isGTF);
}

/*
    destroys the GeneAnnoParser
*/
GeneAnnoParser::~GeneAnnoParser()
{
    // destroy the associated GFFParser
    delete(this->gffParser);
}

/*
    parses through the entire file
*/
GFFData
GeneAnnoParser::parse()
{
    // work out which parsing function should be used
    void (GeneAnnoParser::*parseFunction)(GFFRecord*) =
        this->selectParseFunction();

    // start parsing through
    GFFRecord rec = gffParser->parseNextRecord();
    while (!(gffParser->isEmpty())) {
        if (!rec.broken) {
            (this->*parseFunction)(&rec);
        }
        rec = gffParser->parseNextRecord();
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
        this->annotationSource = gffParser->guessAnnotationSource();
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
    // first check that the record has a gene_id
    if (!rec->hasAttribute("gene_id")) {
            std::cout << "Record did not have 'gene_id' attribute: " << rec->printAttributes() << "\n";
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
    std::cout << "running parseGENCODE, todo\n";
}

/*
    wrapper for parsing a GFF or a GTF into a GFFData object
*/
GFFData
parseGeneAnno(std::string filename)
{
    GeneAnnoParser * geneAnnoParser = new GeneAnnoParser(filename);
    GFFData gffData = geneAnnoParser->parse();
    delete(geneAnnoParser);
    return gffData;
}

/*
    checks whether a filename contains "gtf"
    I guess it could be fooled if your GFF has ".gtf" somewhere in its name...
    but otherwise it'll do
*/
bool
isFileGTF(std::string filename)
{
    return filename.find(".gtf") != std::string::npos;
}