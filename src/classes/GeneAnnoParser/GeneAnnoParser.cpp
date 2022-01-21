#include "GeneAnnoParser.h"

#include <algorithm>
#include <functional>
#include <string>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <iostream>
#include <Rcpp.h>

#include "GFFRecord.h"
#include "GFFParser.h"

#include "../Parser.h"
#include "../Pos.h"
#include "../StartEndPair.h"
#include "../GFFData.h"
// int
// main()
// {
//     GeneAnnoParser * 
//     geneAnnoParser = new GeneAnnoParser("data/isoform_annotated.gff3");
//     geneAnnoParser->parse();
// }

/*****************************************************************************/

GeneAnnoParser::GeneAnnoParser(std::string filename, bool isGTF)
{
    this->filename = filename;
    this->isGTF = isGTF;
    this->annotationSource = "Ensembl"; // this should not be a default, but needs to be decided based on old guess_annotation_source function

}

/*
    parses through the entire file
*/
GFFData
GeneAnnoParser::parse()
{
    // set up the file parser
    GFFParser fileParser(filename, isGTF ? "GTF" : "GFF");

    // work out which parsing function should be used
    // void (GeneAnnoParser::*parseFunction)(GFFRecord*) =
   	ParseFunction parseFunction = selectParseFunction();

    // start parsing through
    GFFRecord rec = fileParser.parseNextRecord();
    while (!fileParser.isEmpty()) {
        if (!rec.broken) {
            // (this->*parseFunction)(&rec);
			parseFunction(&rec);
		}
        rec = fileParser.parseNextRecord();
    }

    this->gffData.removeTranscriptDuplicates(this->isGTF);

	fileParser.close();
    return this->gffData;
}

/*
    selects the appropriate function to point to for parsing,
    based on whether the file is a GTF and what its annotation source is
*/
// ParseFunction
ParseFunction
GeneAnnoParser::selectParseFunction()
{
    // ParseFunction parseFunction;
    ParseFunction parseFunction;
	if (this->isGTF) {
        parseFunction = [this](GFFRecord *rec){this->parseGTF(rec);};
		// parseFunction = &GeneAnnoParser::parseGTF;
    } else {
        this->annotationSource = gffParser->guessAnnotationSource();
        if (this->annotationSource == "Ensembl") {
            parseFunction = [this](GFFRecord *rec){this->parseEnsembl(rec);};
			// parseFunction = &GeneAnnoParser::parseEnsembl;
        } else {
			parseFunction = [this](GFFRecord *rec){this->parseGENCODE(rec);};
            // parseFunction = &GeneAnnoParser::parseGENCODE;
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
            Rcpp::Rcout << "Record did not have 'gene_id' attribute: " << rec->printAttributes() << "\n";
            return;
    }
    std::string gene_id = rec->attributes["gene_id"];

	std::vector<std::string> chr_vec = this->gffData.chr_to_gene[rec->seqname];
    if (std::find(chr_vec.begin(), chr_vec.end(), gene_id) == chr_vec.end()) {
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
	if (rec->hasAttribute("gene_id")) {
		gffData.chr_to_gene[rec->seqname].push_back(rec->attributes["gene_id"]);
	}

	ParseResult pr = parseKeyValue(rec->attributes["Parent"], ':');
	std::string parent_type 	= pr.first;
	std::string parent_gene_id 	= pr.second;
	std::string transcript_id 	= rec->attributes["transcript_id"];

	if (rec->hasAttribute("Parent") && parent_type == "gene") {
		gffData.gene_to_transcript[parent_gene_id].push_back(transcript_id);
		gffData.transcript_dict[transcript_id] =  Pos {rec->seqname, rec->start-1, rec->end, rec->strand, parent_gene_id};
	} else if (rec->feature == "exon") {
		if (parent_type != "transcript") {
			Rcpp::Rcout << "Format Error. What to do?\n";
			return;
		}
		gffData.transcript_to_exon[parent_gene_id].push_back(StartEndPair {rec->start-1, rec->end});
	}

	// old, untested
    // if (rec->feature == "gene") {
    //     // found a gene
    //     this->gffData.chr_to_gene[rec->seqname].push_back(rec->attributes["gene_id"]);
    // } else if (rec->feature == "transcript") {
    //     // found a transcript
    // } else if (rec->feature == "exon") {
    //     // found an exon
    // }
    // auto gene_id = rec->attributes["gene_id"];
    // auto transcript_id = rec->attributes["transcript_id"];

    // if (rec->feature == "transcript") {
    //     auto parent = rec->attributes["Parent"];

    //     this->gffData.gene_to_transcript[gene_id].push_back(transcript_id);
    // }
}

/*
    parses a single record of a GENCODE GFF file
*/
void
GeneAnnoParser::parseGENCODE(GFFRecord * rec)
{
    std::cout << "running parseGENCODE, todo\n";
}