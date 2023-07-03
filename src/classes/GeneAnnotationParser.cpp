#include "GeneAnnotationParser.h"

#include <unordered_map>
#include <string>
#include <vector>
#include <set>
#include "zlib.h"

#include "GFFData.h"
#include "GFFRecord.h"  
#include "types.h"
#include "StartEndPair.h"
#include "Pos.h"
#include "../utility/utility.h"

enum AnnotationSource { gencode, ensembl };
AnnotationSource guessAnnotationSource(const std::string &filename) 
{
    gzFile annotation = gzopen(filename.c_str(), "r");

    const int MAXLEN = 256;
    char buf[MAXLEN];
    AnnotationSource returnSource = ensembl;
    // Only look through the first 1000 lines
    for (int i = 0; i < 1000; i++) {
        if (gzgets(annotation, buf, MAXLEN) == NULL) return ensembl;
        std::string line(buf);
        if (line.find("GENCODE") != std::string::npos) {
            returnSource = gencode;
            break;
        } else if (line.find("1\tEnsembl") != std::string::npos) {
            returnSource = ensembl;
            break;
        }
    }

    gzclose(annotation);
    return returnSource;
}

void
removeTranscriptDuplicates(GFFData * const data, bool updateTranscriptDict=true)
{
    for (auto &[transcript, exons] : data->transcript_to_exon) {
        ranges::sort<StartEndPair>(exons);
        if (exons.size() > 1 && exons[0].start == exons[1].start) {
            std::vector<exon> new_exons = {exons[0]};
            for (const auto &ex : exons) {
                if (ex.start != new_exons.back().start && ex.end != new_exons.back().start)
                    new_exons.push_back(ex);
            }

            data->transcript_to_exon[transcript] = new_exons;
        }

        if (updateTranscriptDict) {
            data->transcript_dict[transcript] = 
                Pos(data->transcript_dict[transcript].chr, data->transcript_to_exon[transcript][0].start,
                    data->transcript_to_exon[transcript].back().end, data->transcript_dict[transcript].strand,
                    data->transcript_dict[transcript].parent_id);
        }
    }

    for (const auto &[gene, transcripts] : data->gene_to_transcript) {
        // remove duplicates from transcript
        // sort transcripts
        std::set<std::string> transcript_set(transcripts.begin(), transcripts.end());
        data->gene_to_transcript[gene] = std::vector<std::string> (transcript_set.begin(), transcript_set.end());
    }
}

GFFData parse_gtf(const std::string &filename)
{
    GFFData data;
    auto attributeParsingFunction = GFFRecord::chooseAttributesFunc(filename);

    gzFile annotationFile = gzopen(filename.c_str(), "r");
    const int MAXLEN = 65535; // maximum char array size
    char buf[MAXLEN];
    while(gzgets(annotationFile, buf, MAXLEN) != NULL) {
        GFFRecord rec = GFFRecord::parseGFFRecord(std::string(buf), attributeParsingFunction);
        
        if (rec.feature == "gene") {
            data.chr_to_gene[rec.seqname].push_back(rec.attributes["gene_id"]);
        } else if (rec.feature == "transcript") {
            if (rec.attributes.count("gene_id") == 0) continue;

            std::string gene_id = rec.attributes["gene_id"];
            if (ranges::doesNotContain(data.chr_to_gene[rec.seqname], gene_id))
                data.chr_to_gene[rec.seqname].push_back(gene_id);
            data.gene_to_transcript[gene_id].push_back(rec.attributes["transcript_id"]);
            data.transcript_dict[rec.attributes["transcript_id"]] = 
                Pos(rec.seqname, rec.start - 1, rec.end, rec.strand, gene_id);
        } else if (rec.feature == "exon") {
            if (rec.attributes.count("gene_id") == 0) continue;

            if (data.chr_to_gene.count(rec.seqname) == 0 
                    || ranges::doesNotContain(data.chr_to_gene[rec.seqname], rec.attributes["gene_id"])) {
                data.chr_to_gene[rec.seqname].push_back(rec.attributes["gene_id"]);
            }
            if (rec.attributes.count("transcript_id")) {
                if (data.transcript_dict.count(rec.attributes["transcript_id"]) == 0) {
                    data.gene_to_transcript[rec.attributes["gene_id"]].push_back(rec.attributes["transcript_id"]);
                    data.transcript_dict[rec.attributes["transcript_id"]] = 
                        Pos(rec.seqname, -1, 1, rec.strand, rec.attributes["gene_id"]);
                }
                data.transcript_to_exon[rec.attributes["transcript_id"]].push_back(
                    StartEndPair(rec.start - 1, rec.end)
                );
            }
        }
        
    }
    gzclose(annotationFile);

    removeTranscriptDuplicates(&data, true);

    return data;
}

GFFData parse_gff(const std::string &filename)
{
    GFFData data;

    AnnotationSource annotation_source = guessAnnotationSource(filename);
    auto attributeParsingFunction = GFFRecord::chooseAttributesFunc(filename);
    
    gzFile annotationFile = gzopen(filename.c_str(), "r");
    const int MAXLEN = 1000;
    char buf[MAXLEN];
    while(gzgets(annotationFile, buf, MAXLEN) != NULL) {
        GFFRecord rec = GFFRecord::parseGFFRecord(std::string(buf), attributeParsingFunction);
        if (annotation_source == ensembl){
            if (rec.attributes.count("gene_id")) 
                data.chr_to_gene[rec.seqname].push_back(rec.attributes["gene_id"]);
            if (rec.attributes.count("Parent") 
                    && rec.attributes["Parent"].find("gene:") != std::string::npos) {
                std::string gene_id = rec.attributes["Parent"].substr(
                    rec.attributes["Parent"].find_first_of(':') + 1);
                data.gene_to_transcript[gene_id].push_back(rec.attributes["transcript_id"]);
                data.transcript_dict[rec.attributes["transcript_id"]] = 
                    Pos (rec.seqname, rec.start - 1, rec.end, rec.strand, gene_id);
            } else if (rec.feature == "exon") {
                if (rec.attributes["Parent"].find("transcript:") == std::string::npos)
                    throw std::invalid_argument("format error");
                std::string gene_id = rec.attributes["Parent"].substr(
                    rec.attributes["Parent"].find_first_of(':') + 1);
                data.transcript_to_exon[gene_id].push_back( StartEndPair(rec.start - 1, rec.end) );
            }
        } else if (annotation_source == gencode) {
            if (rec.feature == "gene") {
                data.chr_to_gene[rec.seqname].push_back(rec.attributes["gene_id"]);
            } else if (rec.feature == "transcript") {
                std::string gene_id = rec.attributes["Parent"];
                if (ranges::doesNotContain(data.chr_to_gene[rec.seqname], gene_id)) 
                    data.chr_to_gene[rec.seqname].push_back(gene_id);
                
                data.gene_to_transcript["gene_id"].push_back(rec.attributes["transcript_id"]);
                data.transcript_dict[rec.attributes["transcript_id"]] = 
                    Pos( rec.seqname, rec.start - 1, rec.end, rec.strand, gene_id);
            } else if (rec.feature == "exon") {
                data.transcript_to_exon[rec.attributes["Parent"]].push_back(
                    StartEndPair(rec.start - 1, rec.end)
                );
            }
        }
    }
    gzclose(annotationFile);

    removeTranscriptDuplicates(&data, false);

    return data;
}


GFFData parse_gff_file(const std::string &gff_filename)
{
    if (gff_filename.find(".gtf") != std::string::npos) {
        return parse_gtf(gff_filename);
    } 

    return parse_gff(gff_filename);
}
