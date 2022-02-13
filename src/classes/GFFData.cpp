#include "GFFData.h"

/*
    iterates through GFFData and removes 
    any duplicate transcripts that may have been added
*/
void
GFFData::removeTranscriptDuplicates(bool updateTranscriptDict) {
    // Remove duplicates from the transcript_to_exon maps
    for (auto tr : this->transcript_to_exon) {
        // it->first is std::string key
        // it->second is std::vector<StartEndPair> list of pairs
        std::sort(tr.second.begin(), tr.second.end(), StartEndPairCompare);
        // if the forward list has at least two elements and start of first element equals start of second
        if (tr.second.size() >= 2 &&
            (tr.second.begin()->start) == (tr.second.begin() + 1)->start) {
            
            std::vector<StartEndPair> new_ex = {*tr.second.begin()};
			// remove duplicates from the old list of exons stored in transcript_to_exon[x]
            for (auto ex : tr.second) {
				StartEndPair last_ex = *(new_ex.end() - 1);
				if (ex.start != last_ex.start && ex.end != last_ex.start) {
					new_ex.push_back(ex);
				}
            }

            this->transcript_to_exon[tr.first] = new_ex;
        }

		if (updateTranscriptDict) {
                // update transcript_dict[this transcript] with a new Pos object of
                // the correct start and end positions of this exon
                this->transcript_dict[tr.first] = Pos {
                    transcript_dict[tr.first].chr,
                    tr.second.begin()->start,
                    (tr.second.end() - 1)->end,
                    transcript_dict[tr.first].strand,
                    transcript_dict[tr.first].parent_id
                };
            }
    }

	// remove duplicates from gene_to_transcript
    for (auto ge : this->gene_to_transcript) {
		std::unordered_map<std::string, int> set;
		std::vector<std::string> new_genes;
		for (auto el : ge.second) {
			set[el] = 1;
		}

		for (auto unique : set) {
			new_genes.push_back(unique.first);
		}
		this->gene_to_transcript[ge.first] = new_genes;
	}
}


/*
    records all the details of GFFData in a file
*/
void
GFFData::log
( 
    std::string filename
)
{
    std::ofstream
    file (filename);

    file << "chr_to_gene: (size " << this->chr_to_gene.size() << ")\n";
    for (const auto & [chr, gene] : this->chr_to_gene) {
        file << "\t" << chr << ": (size " << gene.size() << ") [";
        for (const auto & g : gene) {
            file << g << ", ";
        }
        file << "]\n";
    }

    file << "\ntranscript_dict: (size " << this->transcript_dict.size() << ")\n";
    for (const auto & [tr, pos] : this->transcript_dict) {
        file << "\t" << tr << ":(" 
            << pos.chr << "," 
            << pos.start << "," 
            << pos.end << "," 
            << pos.parent_id << "," 
            << pos.strand << ")\n";
    }

    file << "\ngene_to_transcript: (size " << this->gene_to_transcript.size() << ")\n";
    for (const auto & [gene, transcript] : this->gene_to_transcript) {
        file << "\t" << gene << ": (size " << transcript.size() << ") [";
        for (const auto & tr : transcript) {
            file << tr << ", ";
        }
        file << "]\n";
    }

    file << "\ntranscript_to_exon: (size " << this->transcript_to_exon.size() << ")\n";
    for (const auto & [transcript, exon] : this->transcript_to_exon) {
        file << "\t" << transcript << ": (size " << exon.size() << ") [";
        for (const auto & ex : exon) {
            file << "(" << ex.start << "," << ex.end << "), ";
        }
        file << "]\n";
    }
}