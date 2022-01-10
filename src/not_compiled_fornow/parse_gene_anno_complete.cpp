#include "parse_gene_anno_complete.hpp"

bool
get_is_GTF(std::string filename)
{
    return filename.find(".gtf") != std::string::npos;
}

std::string
get_annotation_source(std::string filename)
{
    /*
        parse the file and check for "GENCODE" or "Ensembl"
    */
    int line_num = 0;

    std::ifstream file (filename);
    std::string line;
    while (std::getline(file, line)) {
        if (line.find("GENCODE") != std::string::npos) {
            file.close();
            return "GENCODE";
        } else if (line.find("1\tEnsembl") != std::string::npos) {
            file.close();
            return "Ensembl";
        }

        // only check the first 1000 lines
        if (line_num++ > 1000) {
            break;
        }
    }

    // return "Ensembl" by default
    file.close();
    return "Ensembl";
}

static bool vectorContains(std::string a, std::vector<std::string> v) {
	for (auto it : v) {
		if (a == it) {
			return true;
		}
	}
	return false;
}

GFFData
parse_GFF_or_GTF(std::string filename)
{
    GFFData gff_data;

    // check if the file is GTF
    bool is_GTF = get_is_GTF(filename);
    // create a new parser and parse the file
    Parser * parser = new Parser(filename, is_GTF);

    return is_GTF ? 
        parse_GTF(&(parser->records)) : 
        parse_GFF(&(parser->records), get_annotation_source(filename));
}

GFFData
parse_GTF(std::vector<Record> * records)
{
    GFFData gff_data;

    for (auto & rec : *records) {
        // first check that the record has a gene_id
        if (!rec.has_attribute("gene_id")) {
                Rcpp::warning("Record did not have 'gene_id' attribute: %s", rec.format_attributes());
                continue;
        }
        auto gene_id = rec.attributes["gene_id"];
        if (!vectorContains(gene_id, gff_data.chr_to_gene[rec.seqname])) {
            gff_data.chr_to_gene[rec.seqname].push_back(gene_id);
        }


        if (rec.feature == "gene") {
            // found a gene
            // skip the rest of the work
            continue;
        }

        if (!rec.has_attribute("transcript_id")) {
            Rcpp::warning("Feature did not have 'transcript_id' attribute: %s", rec.format_attributes());
        }    
        auto transcript_id = rec.attributes["transcript_id"];
        gff_data.gene_to_transcript[gene_id].push_back(transcript_id);
        StartEndPair start_end = rec.feature == "transcript" ? 
            (StartEndPair){rec.start-1, rec.end} : 
            (StartEndPair){-1, 1};
        gff_data.transcript_dict[transcript_id] = {
            rec.seqname,
            rec.start - 1,
            rec.end,
            rec.strand,
            gene_id
        };

        if (rec.feature == "transcript") {
            // found a transcript
            // no more processing to do
            continue;
        }

        // found an exon
        gff_data.transcript_to_exon[transcript_id].push_back({
            rec.start - 1,
            rec.end
        });
    }
    return gff_data;
}

GFFData
parse_GFF(std::vector<Record> * records, std::string source)
{
    GFFData gff_data;

    if (source == "Ensembl") {
        for (auto & rec : *records) {
            
        }
    }
}

void
parse_one_Ensembl(GFFData * gff_data, Record rec)
{
    if (rec.feature == "gene") {
        // found a gene
        (*gff_data).chr_to_gene[rec.seqname].push_back(rec.attributes["gene_id"]);
    } else if (rec.feature == "transcript") {
        // found a transcript
    } else if (rec.feature == "exon") {
        // found an exon
    }
    auto gene_id = rec.attributes["gene_id"];
    auto transcript_id = rec.attributes["transcript_id"];

    if (rec.feature == "transcript") {
        auto parent = rec.attributes["Parent"];

        (*gff_data).gene_to_transcript[gene_id].push_back(transcript_id);
    }
}

/*****************************************************************************/

/*
    exports the attributes as a formatted string
*/
std::string
Record::format_attributes()
{
    std::string
    output;
    for (const auto & [key, val] : this->attributes) {
        output += key + "\"" + val + "\"" + ";";
    }
    return output;
}