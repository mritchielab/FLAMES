#include "gtf_to_bed.h"

// [[Rcpp::export]]
void
gtf_to_bed_cpp(std::string in_gtf, std::string out_bed, std::string chrom_sizes_file)
{
  std::cout << "started gtf_to_bed_cpp\n";

  std::ifstream gtf (in_gtf);
  std::ofstream bed (out_bed);
  std::ifstream chrom_sizes;

  std::cout << "opened the files\n";

  bool is_bed = (out_bed.substr(out_bed.length() - 3, out_bed.length()) != "psl") ? true : false;

  if (chrom_sizes_file != "") {
    chrom_sizes.open(chrom_sizes_file); 
  }

  std::map<std::string, int>
  chrom_to_size = {};

  std::string line;
  if (chrom_sizes_file != "") {
    while (std::getline(chrom_sizes, line)) {
      std::istringstream stream (line);
      
      std::vector<std::string> tokens;
      std::string token;
      while (std::getline(stream, token, '\t')) {
        tokens.push_back(token);
      }

      chrom_to_size[tokens[0]] = std::stoi(tokens[1]);
    }
  }

  std::set<std::string>
  missing_chroms;

  // lambda function to convert a vector to a comma-separated stringstream
  auto vector_to_str = [] (
    std::vector<int> vector
  ) 
  {
    std::stringstream stringstream = (std::stringstream)("");
    for (auto i : vector) {
      stringstream << i << ',';
    }
    return stringstream.str();
  };

  // lambda function to write rows to a bed file
  // format described here: https://en.wikipedia.org/wiki/BED_(file_format)
  auto bed_write_row = [] (
    std::ofstream& file, 
    std::string chrom,
    int chromStart,
    int chromEnd,
    std::string name,
    int score,
    std::string strand,
    int thickStart,
    int thickEnd,
    std::string itemRgb,
    int blockCount,
    std::string blockSizes,
    std::string blockStarts
  )
  {
    file << chrom << "\t"
        << chromStart << "\t"
        << chromEnd << "\t"
        << name << "\t"
        << score << "\t"
        << strand << "\t"
        << thickStart << "\t"
        << thickEnd << "\t"
        << itemRgb << "\t"
        << blockCount << "\t"
        << blockSizes << "\t"
        << blockStarts << "\n";
  };

  // lambda function to write rows to a psl file
  // format described here: http://genome.ucsc.edu/FAQ/FAQformat#format2
  auto
  psl_write_row = [&bed] (
    int matches,
    int misMatches,
    int repMatches,
    int nCount,
    int qNumInsert,
    int qBaseInsert,
    int tNumInsert,
    int tBaseInsert,
    std::string strand,
    std::string qName,
    int qSize,
    int qStart,
    int qEnd,
    std::string tName,
    int tSize,
    int tStart,
    int tEnd,
    int blockCount,
    std::string blockSizes,
    std::string qStarts,
    std::string tStarts
  ) 
  {
    bed << matches << "\t"
        << misMatches << "\t"
        << repMatches << "\t"
        << nCount << "\t"
        << qNumInsert << "\t"
        << qBaseInsert << "\t"
        << tNumInsert << "\t"
        << qBaseInsert << "\t"
        << tNumInsert << "\t"
        << tBaseInsert << "\t"
        << strand << "\t"
        << qName << "\t"
        << qSize << "\t"
        << qStart << "\t"
        << qEnd << "\t"
        << tName << "\t"
        << tSize << "\t"
        << tStart << "\t"
        << tEnd << "\t"
        << blockCount << "\t"
        << blockSizes << "\t"
        << qStarts << "\t"
        << tStarts << "\n";
  };
  
  // generic function to write any amount of data to a tab-separated CSV row
  auto
  csv_write_row = [] (std::ofstream& file, std::vector<std::string> data)
  {
    for (auto entry : std::vector<std::string>(data.begin(), data.end() - 1)) {
      file << entry << '\t';
    }
    file << data.back() << '\n';
  };

  std::string prev_transcript = "";
  std::string prev_chrom;
  std::string prev_gene;
  std::string prev_strand;

  int blockcount;
  std::vector<int> blockstarts = {};
  std::vector<int> blocksizes = {};

  std::string gtf_line;
  std::string this_transcript;
  int tstart;
  int tend;

  std::string chrom;
  std::string strand;
  int start;
  int end;

  std::vector<int> qstarts;

  while (std::getline(gtf, gtf_line)) {
    if (gtf_line[0] == '#') {
      continue;
    }
    
    std::stringstream gtf_line_stream (gtf_line);
    std::vector<std::string> values;
    std::string value;
    while (std::getline(gtf_line_stream, value, '\t')) {
      values.push_back(value);
    }

    auto chrom = values[0];
    auto ty = values[2];
    auto start = std::stoi(values[3]) - 1;
    auto end = std::stoi(values[4]);
    auto strand = values[6];

    if (ty != "exon") {
      continue;
    }
    // just get a substring of the attributes entry
    auto transcript_id = values[8].find("transcript_id") + 15;
    this_transcript = values[8].substr(transcript_id, values[8].length());
    this_transcript = this_transcript.substr(0, this_transcript.find("\""));
    // once all the exons for a transcript are read, write the psl/bed entry
    if (this_transcript != prev_transcript) {
      if (prev_transcript != "") {
        blockcount = blockstarts.size();
        if (blockcount > 1 && blockstarts[0] > blockstarts[1]) { // we need to reverse the exons
          std::reverse(blocksizes.begin(), blocksizes.end());
          std::reverse(blockstarts.begin(), blockstarts.end());
        }
        // target (eg chrom)
        tstart = blockstarts.front();
        tend = blockstarts.back() + blockstarts.back();

        // query (eg transcript)
        auto qsize = 0;
        for (auto num : blocksizes) {
          qsize += num;
        }

        auto qname = prev_transcript; 
        qstarts = {};

        if (!(is_bed)) { // psl specific
          auto pos = 0;
          qstarts = {pos};

          for (auto b : std::vector<int>(blocksizes.begin(), blocksizes.end() - 1)) {
            pos += b;
            qstarts.push_back(pos);
          }
        }
        if (is_bed) { // bed specific
          std::vector<int>
          relative_blockstarts = {}; // stores the block starts relative to the transcript start
          for (auto block : blockstarts) {
            relative_blockstarts.push_back(block - tstart);
          }

          // log with bed formatting
          bed <<prev_chrom<<"\t"
                    <<tstart<<"\t"
                    <<tend<<"\t"
                    <<qname<<"\t"
                    <<1000<<"\t"
                    <<prev_strand<<"\t"
                    <<tstart<<"\t"
                    <<tend<<"\t"
                    <<0<<"\t"
                    <<blockcount<<"\t"
                    <<vector_to_str(blocksizes)<<"\t"
                    <<vector_to_str(relative_blockstarts)<<"\n";
          
          // bed_write_row(bed, prev_chrom, tstart, tend, qname, 1000, prev_strand, tstart, tend, 0, blockcount, vector_to_str(blocksizes), vector_to_str(relative_blockstarts));
        } else { // psl specific

          // to get the query sequence size value for the pcl, we look in the dictionary
          int tsize = 0;
          if (chrom_to_size.size() > 0) {
            if (chrom_to_size.count(prev_chrom) > 0) {
              tsize = chrom_to_size[prev_chrom];
            } else {
              missing_chroms.insert(prev_chrom);
            }
          }

          psl_write_row(0, 0, 0, 0, 0, 0, 0, 0, prev_strand, qname, qsize, 0, qsize, prev_chrom, tsize, tstart, tend, (int)blockcount, vector_to_str(blocksizes), vector_to_str(qstarts), vector_to_str(blockstarts));
        }
      }

      // reset everything for the next line
      blockstarts = {};
      blocksizes = {};
      prev_transcript = this_transcript;

      // just get the gene_id from the attributes entry and assign it to prev_gene
      auto gene_id = values[8].find("gene_id") + 9;
      prev_gene = values[8].substr(gene_id, values[8].length());
      prev_gene = prev_gene.substr(0, this_transcript.find("\""));

      prev_chrom = chrom;
      prev_strand = strand;
    }

    blockstarts.push_back(start);
    blocksizes.push_back(end - start);
  }
  // last entry...
  if (blockcount > 1 && blockstarts[0] > blockstarts[1]) { // need to reverse exons
    std::reverse(blocksizes.begin(), blocksizes.end());
    std::reverse(blockstarts.begin(), blockstarts.end());
  }

  // target (eg chrom)
  tstart = blockstarts.front();
  tend = blockstarts.back() + blockstarts.back();

  // query (eg transcript)
  auto qsize = 0;
  for (auto num : blocksizes) {
    qsize += num;
  }

  auto qname = this_transcript;

  if (is_bed) { // bed specific
    std::vector<int>
    relative_blockstarts; // stores the block starts relative to the transcript start
    for (auto block : blockstarts) {
      relative_blockstarts.push_back(block - tstart);
    }

    // bed_write_row(bed, chrom, tstart, tend, qname, 1000, strand, start, end, 0, blockcount, vector_to_str(blocksizes), vector_to_str(relative_blockstarts));
    // log with bed formatting
    bed <<chrom<<"\t"
              <<tstart<<"\t"
              <<tend<<"\t"
              <<qname<<"\t"
              <<1000<<"\t"
              <<strand<<"\t"
              <<start<<"\t"
              <<end<<"\t"
              <<0<<"\t"
              <<blockcount<<"\t"
              <<vector_to_str(blocksizes)<<"\t"
              <<vector_to_str(relative_blockstarts)<<"\n";

  } else { // psl specific
    // to get the query sequence size value for the pcl, we look in the dictionary
    int tsize = 0;
    if (chrom_to_size.size() > 0) {
      if (chrom_to_size.count(prev_chrom) > 0) {
        tsize = chrom_to_size[prev_chrom];
      } else {
        missing_chroms.insert(prev_chrom);
      }
    }

    psl_write_row(0, 0, 0, 0, 0, 0, 0, 0, prev_strand, qname, qsize, 0, qsize, prev_chrom, tsize, tstart, tend, int(blockcount), vector_to_str(blocksizes), vector_to_str(qstarts), vector_to_str(blockstarts));
  }

  if (missing_chroms.size() > 0) {
    std::cout << "chromosomes found in gtf but not in chrom_sizes file: ";
    for (auto chrom : missing_chroms) {
      std::cout << chrom;
    }
    std::cout << "\n";
  }
  std::cout << "finished gtf_to_bed_cpp\n";
}