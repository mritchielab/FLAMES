// Copyright 2022 Nadia Davidson 
// This program is distributed under the MIT License.
// We also ask that you cite this software in publications
// where you made use of it for any part of the data analysis.

#include <iostream>
#include <fstream>
#include <istream>
#include <string>
#include <sstream>
#include <unordered_set>
#include <unordered_map>
#include <stdlib.h>
#include <unistd.h>
#include <vector>
#include <algorithm>
#include <thread>
#include <numeric>
#include <Rcpp.h>
#include <ctime>
#include <filesystem>
// [[Rcpp::plugins(cpp17)]]

#include "./utility/edlib-1.2.7/edlib.h"


const static std::string VERSION="0.96.2";


// compliment nucleotides - used to reverse compliment string
char compliment(char& c){
  switch(c){
  case 'A' : return 'T';
  case 'T' : return 'A';
  case 'G' : return 'C';
  case 'C' : return 'G';
  default: return 'N';
  }
}

// Reverse complement
std::string reverse_compliment(const std::string& seq) {
  std::string rev_seq(seq.rbegin(), seq.rend());
  std::transform(rev_seq.begin(), rev_seq.end(), rev_seq.begin(),
                 [](char c) { return compliment(c); });
  return rev_seq;
}

//Holds the search string patterns
struct SearchSeq {
  std::string primer;
  std::string polyA;
  std::string umi_seq;
  std::string temp_barcode;
} search_pattern ;

//Holds the found barcode and associated information 
struct Barcode {
  std::string barcode;
  std::string umi;
  int editd;
  int flank_editd;
  int flank_start;
  int flank_end;
  bool unambiguous;
} ;

struct SearchResult {
  std::string read_id;
  std::string qual_scores;
  std::string line;
  std::string rev_line;
  std::vector<Barcode> vec_bc_for;
  std::vector<Barcode> vec_bc_rev;
};

// Code for fast edit distance calculation for short sequences modified from
// https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance#C++
// s2 is always assumned to be the shorter string (barcode)
unsigned int edit_distance(const std::string& s1, const std::string& s2, unsigned int &end, int max_editd){

  const std::string_view s1_view(s1);
  const std::string_view s2_view(s2);

  const std::size_t len1 = s1_view.size() + 1;
  const std::size_t len2 = s2_view.size() + 1;

  std::vector<unsigned int> dist_holder(len1 * len2);
  //initialise the edit distance matrix.
  //penalise for gaps at the start and end of the shorter sequence (j)
  //but not for shifting the start/end of the longer sequence (i,0)
  dist_holder[0] = 0; //[0][0]
  for (std::size_t j = 1; j < len2; ++j)
    dist_holder[j] = j; //[0][j];
  for (std::size_t i = 1; i < len1; ++i)
    dist_holder[i * len2] = 0; //[i][0];

  unsigned int best = len2;
  end = len1 - 1;

  // loop over the distance matrix elements and calculate running distance
  for (std::size_t j = 1; j < len2; ++j) {
    bool any_below_threshold = false; // flag used for early exit
    for (std::size_t i = 1; i < len1; ++i) {
      unsigned int sub = (s1_view[i - 1] == s2_view[j - 1]) ? 0 : 1; // match / mismatch score

      const unsigned int &top_left = dist_holder[(i - 1) * len2 + (j - 1)];
      const unsigned int &left = dist_holder[i * len2 + (j - 1)];
      const unsigned int &top = dist_holder[(i - 1) * len2 + j];

      unsigned int min_value = std::min({top + 1, left + 1, top_left + sub});
      dist_holder[i * len2 + j] = min_value;

      if (min_value <= max_editd)
        any_below_threshold = true;
      if (j == (len2 - 1) && min_value < best) {
        // if this is the last row in j
        // check if this is the best running score
        // update the end position of alignment
        best = min_value;
        end = i;
      }
    }
    if(!any_below_threshold){ //early exit to save time.
      return(100);
    }
  }
  return best; //return edit distance
}


//Given a string 'seq' search for substring with primer and polyT sequence followed by
//a targeted search in the region for barcode
//Seaquence seearch is performed using edlib

Barcode get_barcode(const std::string & seq,
		    const std::unordered_set<std::string> &known_barcodes,
		    int flank_max_editd,
		    int barcode_max_editd){
  
  const int OFFSET=5; //wiggle room in bases of the expected barcode start site to search.

  //initialise struct variables for return:
  Barcode barcode;
  barcode.editd=100; barcode.flank_editd=100; barcode.unambiguous = false;

  //initialise edlib configuration
  EdlibEqualityPair additionalEqualities[5] = {{'?','A'},{'?','C'},{'?','G'},{'?','T'},{'?','N'}};
  EdlibAlignConfig edlibConf = {flank_max_editd, EDLIB_MODE_HW, EDLIB_TASK_PATH, additionalEqualities, 5};

  //search for primer and ployT (barcode and umi as wildcards)
  std::string search_string=
    search_pattern.primer+
    search_pattern.temp_barcode+
    search_pattern.umi_seq+
    search_pattern.polyA;
  EdlibAlignResult result = edlibAlign(search_string.c_str(), search_string.length(), seq.c_str(), seq.length(), edlibConf);
  if(result.status != EDLIB_STATUS_OK || result.numLocations==0 ){
    edlibFreeAlignResult(result);
    return(barcode); // no match found - return
  } //fill in info about found primer and polyT location
  barcode.flank_editd=result.editDistance;
  barcode.flank_start=result.startLocations[0];
  barcode.flank_end=seq.find_first_not_of('T', result.endLocations[0]);

  // Extract sub-patterns from aligment directly
  std::vector<long unsigned int> subpattern_lengths = {
    search_pattern.primer.length(),
    search_pattern.temp_barcode.length(),
    search_pattern.umi_seq.length(),
    search_pattern.polyA.length()
  };

  std::vector<long unsigned int> subpattern_ends;
  subpattern_ends.resize(subpattern_lengths.size());
  std::inclusive_scan(subpattern_lengths.begin(), subpattern_lengths.end(), subpattern_ends.begin());

  std::vector<int> read_to_subpatterns;
  read_to_subpatterns.reserve(subpattern_ends.size());

  // initialise pointers
  int i_read = barcode.flank_start;
  int i_pattern = 0;
  int i_subpattern = 0;

  // walk through edlib aligment
  // 0 for match
  // 1 for insertion to target
  // 2 for insertion to query
  // 3 for mismatch 
  std::vector<unsigned char> alignment_vector(result.alignment, result.alignment + result.alignmentLength);
  for (const auto& value : alignment_vector) {
    if (value != 1) {
      i_read++;
    }
    if (value != 2) {
      i_pattern++;
    }
    if (i_pattern >= subpattern_ends[i_subpattern]) {
      read_to_subpatterns.emplace_back(i_read);
      i_subpattern++;
    }
  }

  edlibFreeAlignResult(result);
  
  //if not checking against known list of barcodes, return sequence after the primer
  //also check for a perfect match straight up as this will save computer later.
  std::string exact_bc=seq.substr(read_to_subpatterns[0], read_to_subpatterns[1] - read_to_subpatterns[0]);
  if(known_barcodes.empty() || (known_barcodes.find(exact_bc) != known_barcodes.end())){ 
    barcode.barcode=exact_bc;
    barcode.editd=0;
    barcode.unambiguous=true;
    barcode.umi=seq.substr(read_to_subpatterns[1], search_pattern.umi_seq.length());
    return(barcode);
  }
  
  // otherwise widen our search space and the look for matches with errors
  std::string barcode_seq=seq.substr(read_to_subpatterns[0]-OFFSET,search_pattern.temp_barcode.length()+2*OFFSET);
  
  //iterate over all the known barcodes, checking each sequentially
  unsigned int editDistance; unsigned int endDistance;
  for (const auto &known_bc : known_barcodes) {
    editDistance = edit_distance(barcode_seq,known_bc,endDistance,barcode_max_editd);
    if (editDistance == barcode.editd) {
      barcode.unambiguous=false;    
    } else if (editDistance < barcode.editd && editDistance <= barcode_max_editd) { //if best so far, update
      barcode.unambiguous=true;
      barcode.editd=editDistance;
      barcode.barcode=known_bc;
      barcode.umi=seq.substr(read_to_subpatterns[0]-OFFSET+endDistance,search_pattern.umi_seq.length());//assumes no error in UMI seq.
      if(editDistance==0) { //if perfect match is found we're done.
        return(barcode);
      }
    }
  }
  return(barcode); //return the best matched barcode and associated information
}

//search a read for one or more barcodes (parent function that calls get_barcode)
std::vector<Barcode> big_barcode_search(const std::string & sequence,const std::unordered_set<std::string> & known_barcodes,
				   int max_flank_editd, int max_editd){
  std::vector<Barcode> return_vec; //vector of all the barcodes found

  //search for barcode
  Barcode result=get_barcode(sequence,known_barcodes,max_flank_editd,max_editd); //,ss);
  if(result.editd<=max_editd && result.unambiguous) //add to return vector if edit distance small enough
    return_vec.emplace_back(result);
  
  //if a result was found, mask out the flanking sequence and search again in case there are more.
  if (!return_vec.empty()) {
    std::string masked_sequence = sequence;
    for(const auto &barcode : return_vec){
      int flank_length;
      if (barcode.flank_end == std::string::npos) {
        flank_length = masked_sequence.length() - 1 - barcode.flank_start;
      } else {
        flank_length = barcode.flank_end - barcode.flank_start;
      }
      masked_sequence.replace(barcode.flank_start, flank_length,std::string(flank_length,'X'));
    } //recursively call this function until no more barcodes are found
  std::vector<Barcode> masked_res = big_barcode_search(masked_sequence,known_barcodes,max_flank_editd,max_editd);
    return_vec.insert(return_vec.end(),masked_res.begin(),masked_res.end()); //add to result
  }
    
  return(return_vec);
}

// print information about barcodes
void print_stats(const std::string& read_id, const std::vector<Barcode>& vec_bc, std::ostream& out_stream) {
  for (const auto& barcode : vec_bc) {
    out_stream << read_id << '\t'
               << barcode.barcode << "\t"
	             << barcode.flank_editd << "\t"
	             << barcode.editd << "\t"
	             << barcode.umi << "\t"
 	             << (barcode.flank_end == std::string::npos ? "True" : "False") << "\n";
  }
}

void print_line(const std::string& id, const std::string& read, const std::string& quals, std::ostream& out_stream) {

  const char delimiter = quals.empty() ? '>' : '@';
  
  //output to the new read file
  out_stream << delimiter << id << '\n'
             << read << '\n';
  
  if (!quals.empty()) {
    out_stream << '+' << id << '\n'
               << quals << '\n';
  }
}

//print fastq or fasta lines..
void print_read(const std::string &read_id, const std::string &read,
                const std::string &qual, const std::vector<Barcode> &vec_bc,
                std::ofstream &outstream, 
                std::unordered_set<std::string> &found_barcodes,
                bool trim_barcodes) {
  // loop over the barcodes found... usually will just be one
  for (int b = 0; b < vec_bc.size(); b++) {

    // format the new read id. Using FLAMES format.
    std::ostringstream ss;
    ss << (b + 1) << "of" << vec_bc.size();
    const std::string &barcode = vec_bc[b].barcode;
    std::string new_read_id =
        barcode + "_" + vec_bc[b].umi + "#" + read_id + ss.str();

    // work out the start and end base in case multiple barcodes
    if (vec_bc.at(b).flank_end == std::string::npos) {
      continue;
    }
    int read_start = vec_bc[b].flank_end;
    int read_length = read.length() - read_start;
    for (int f = 0; f < vec_bc.size(); f++) {
      int temp_read_length = vec_bc[f].flank_start - read_start;
      if (temp_read_length > 0 && temp_read_length < read_length)
        read_length = temp_read_length;
    }
    std::string qual_new = (qual.empty() ? std::string() :
                                 qual.substr(read_start, read_length));
    std::string read_new = read.substr(read_start, read_length);

    if (b == 0 && !trim_barcodes) { // override if read shouldn't be cut
      new_read_id = read_id;
      read_new = read;
      qual_new = qual;
      b = vec_bc.size(); // force loop to exit after this iteration
    }

    print_line(new_read_id, read_new, qual_new, outstream);
  }
}

// separated out from main so that this can be run with threads
void search_read(std::vector<SearchResult> & reads, std::unordered_set<std::string> & known_barcodes,
			 int flank_edit_distance, int edit_distance){
  for (auto& read : reads) {
    //forward search
    read.vec_bc_for = big_barcode_search(read.line, known_barcodes, flank_edit_distance, edit_distance);

    //Check the reverse compliment of the read
    read.rev_line = reverse_compliment(read.line);
    read.vec_bc_rev = big_barcode_search(read.rev_line, known_barcodes, flank_edit_distance, edit_distance);
  }
}

//' @export
// [[Rcpp::export]]
int flexiplex(Rcpp::String reads_in, Rcpp::String barcodes_file, bool bc_as_readid, int max_bc_editdistance,
              int max_flank_editdistance, Rcpp::StringVector pattern, Rcpp::String reads_out,
              Rcpp::String stats_out, Rcpp::String bc_out, int n_threads){
  std::ios_base::sync_with_stdio(false);

  bool remove_barcodes=true;
  search_pattern.primer = "CTACACGACGCTCTTCCGATCT"; //(p)
  search_pattern.polyA = std::string(9,'T'); //(T)
  search_pattern.umi_seq = std::string(12,'?'); //(length u)
  search_pattern.temp_barcode = std::string(16,'?'); //(length b)

  for (auto key: std::vector<std::string>{"primer", "polyT", "umi_seq", "barcode_seq"}) {
    if (pattern.containsElementNamed(key.c_str())) {
      if (key == "primer") {
        search_pattern.primer = pattern[key];
      } else if (key == "polyT") {
        search_pattern.polyA = pattern[key];
      } else if (key == "umi_seq") {
        search_pattern.umi_seq = pattern[key];
      } else if (key == "barcode_seq") {
        search_pattern.temp_barcode = pattern[key];
      }
    }
  }
 
  Rcpp::Rcout << "FLEXIPLEX " << VERSION << "\n";
  Rcpp::Rcout << "Setting max barcode edit distance to "<< max_bc_editdistance << "\n";
  Rcpp::Rcout << "Setting max flanking sequence edit distance to "<< max_flank_editdistance << "\n";
  Rcpp::Rcout << "Setting read IDs to be " << (remove_barcodes ? "" : "not") << " replaced" << "\n";
  Rcpp::Rcout << "Setting number of threads to "<< n_threads << "\n";

  //Set of known barcodes 
  std::unordered_set<std::string> known_barcodes;
  std::unordered_set<std::string> found_barcodes;

  /**** READ BARCODE LIST FROM FILE ******/
  if (barcodes_file != "") {
    std::string bc;
    std::ifstream bc_file(barcodes_file);
    if (bc_file.is_open()) {
      Rcpp::Rcout << "Setting known barcodes from "<< barcodes_file.get_cstring() << "\n";
      std::string line;
      while ( getline (bc_file,line) ){
        std::istringstream line_stream(line);
        line_stream >> bc;
        known_barcodes.insert(bc); 
      }
    } else { //if the string given isn't a file
      std::stringstream bc_list(barcodes_file);
	    while (getline(bc_list, bc, ',')){ //tokenize
        known_barcodes.insert(bc);
      }
    }
    bc_file.close();
    Rcpp::Rcout << "Number of known barcodes: " << known_barcodes.size() << "\n";
  }

  std::ifstream ifreads_in(reads_in);
  if(!(ifreads_in.is_open())){
    Rcpp::stop(std::string("Unable to open file ") + reads_in.get_cstring());
  }

  std::ofstream outstream(reads_out, std::ios_base::app);

  /********* FIND BARCODE IN READS ********/
  std::string sequence;
  int bc_count=0;
  int r_count=0;
  int multi_bc_count=0;
  
  std::ofstream out_stat_file;

  if(known_barcodes.size()>0){
    if (std::filesystem::exists(stats_out.get_cstring())) {
      out_stat_file.open(stats_out, std::ios_base::app);
    } else {
      out_stat_file.open(stats_out);
      out_stat_file << "Read\tCellBarcode\tFlankEditDist\tBarcodeEditDist\tUMI\tTooShort"<<"\n";
    }
  }
  Rcpp::Rcout << "Searching for barcodes..." << "\n";
  bool is_fastq=true;
  std::unordered_map< std::string, int > barcode_counts; 
  std::string read_id_line;
  if(getline (ifreads_in,read_id_line)){ //check the first line for file type
    if(read_id_line[0]=='>'){ is_fastq=false;
    } else if (read_id_line[0] == '@'){ //fasta
    } else {
      Rcpp::stop("Unknown read format... exiting");
    }
  }
  
  std::string line;
  while ( getline (ifreads_in,line) ){
    const int buffer_size = 2000; //number of reads to pass to one thread.
    std::vector<std::vector<SearchResult>> sr_v(n_threads);
    for(int i=0; i<n_threads; i++)
      sr_v[i]=std::vector<SearchResult>(buffer_size); 
    std::vector<std::thread> threads(n_threads);
    for(int t=0; t < n_threads; t++){ //get n_threads*buffer number or reads..
      for(int b=0 ; b < buffer_size ; b++){ 
	SearchResult & sr = sr_v[t][b];
	sr.line=line;
  std::string read_id;
	//sr.read_id= read_id_line.substr(1,read_id_line.find_first_not_of(" \t")-1);      
  std::istringstream line_stream(read_id_line);
	line_stream >> sr.read_id;
	sr.read_id.erase(0,1);

	    
	//      string qual_scores="";
	if(!is_fastq){ //fastq (account for multi-lines per read)
    std::string buffer_string; 
	  while(getline(ifreads_in,buffer_string) && buffer_string[0]!='>')
	    sr.line+=buffer_string;
	  read_id_line=buffer_string;
	} else { //fastq (get quality scores)
	  for(int s=0; s<2; s++) getline(ifreads_in,sr.qual_scores);
	  getline(ifreads_in,read_id_line);
	}
      
	r_count++; //progress counter
	if(r_count % 100000 == 0)
	  Rcpp::Rcout << r_count/((double) 1000000 ) << " million reads processed.." << "\n";

	//this is quite ugly, must be a better way to do this..
	if(b==buffer_size-1 && t==n_threads-1){
	  break; //if it's the last in the chunk don't getline as this happens in the while statement
	} else if( !getline(ifreads_in,line)){ //case we are at the end of the reads.
	  sr_v[t].resize(b+1);
	  threads[t]=std::thread(search_read,ref(sr_v[t]),ref(known_barcodes),max_flank_editdistance,max_bc_editdistance);
	  for(int t2=t+1; t2 < n_threads ; t2++) sr_v[t2].resize(0);
	  goto print_result; //advance the line
	}
      }
      // send reads to the thread
      threads[t]=std::thread(search_read,ref(sr_v[t]),ref(known_barcodes),max_flank_editdistance,max_bc_editdistance);
    }
  print_result:
    
    for(int t=0; t < sr_v.size(); t++){ //loop over the threads and print out ther results
      if(sr_v[t].size()>0) threads[t].join(); // wait for the threads to finish before printing
      
      for(int r=0; r< sr_v[t].size(); r++){ // loop over the reads
	
	for(int b=0; b<sr_v[t][r].vec_bc_for.size(); b++)
	  barcode_counts[sr_v[t][r].vec_bc_for.at(b).barcode]++;
	for(int b=0; b<sr_v[t][r].vec_bc_rev.size(); b++)
	  barcode_counts[sr_v[t][r].vec_bc_rev.at(b).barcode]++;
	
	if((sr_v[t][r].vec_bc_for.size()+sr_v[t][r].vec_bc_rev.size())>0)
	  bc_count++;
	if((sr_v[t][r].vec_bc_for.size()+sr_v[t][r].vec_bc_rev.size())>1 ){
	  multi_bc_count++;
	}
	
	if(known_barcodes.size()!=0){ // if we are just looking for all possible barcodes don't output reads etc.
	  
	  print_stats(sr_v[t][r].read_id, sr_v[t][r].vec_bc_for, out_stat_file);
	  print_stats(sr_v[t][r].read_id, sr_v[t][r].vec_bc_rev, out_stat_file);
	  
	  print_read(sr_v[t][r].read_id+"_+",sr_v[t][r].line,sr_v[t][r].qual_scores,sr_v[t][r].vec_bc_for,
		     outstream,found_barcodes,remove_barcodes);
	  reverse(sr_v[t][r].qual_scores.begin(),sr_v[t][r].qual_scores.end());
	  if(remove_barcodes || sr_v[t][r].vec_bc_for.size()==0) //case we just want to print read once if multiple bc found.
	    print_read(sr_v[t][r].read_id+"_-",sr_v[t][r].rev_line,sr_v[t][r].qual_scores,sr_v[t][r].vec_bc_rev,
		       outstream,found_barcodes,remove_barcodes);
	}
      }
    }
  }
  ifreads_in.close();
  
  Rcpp::Rcout << "Number of reads processed: " << r_count << "\n";
  Rcpp::Rcout << "Number of reads where a barcode was found: " << bc_count << "\n";
  Rcpp::Rcout << "Number of reads where more than one barcode was found: " << multi_bc_count << "\n";
  Rcpp::Rcout << "All done!" << "\n";

  if(known_barcodes.size()>0){
    out_stat_file.close();
    return(0);
  }
  
  if(barcode_counts.size()==0)
    return(0);
  
  typedef std::pair<std::string, int> pair;
  std::vector<pair> bc_vec;
  copy(barcode_counts.begin(),barcode_counts.end(), std::back_inserter<std::vector<pair>>(bc_vec));
  sort(bc_vec.begin(), bc_vec.end(),[](const pair &l, const pair &r){
    if (l.second != r.second)
      return l.second > r.second;
    return l.first < r.first;
  });
  std::vector<int> hist(bc_vec[0].second);
  std::ofstream out_bc_file;
  out_bc_file.open(bc_out);
  for (auto const &bc_pair: bc_vec){
    out_bc_file << bc_pair.first << "\t" << bc_pair.second << "\n";
    hist[bc_pair.second-1]++;
  }
  out_bc_file.close();

  Rcpp::Rcout << "Reads\tBarcodes" << "\n";
  for(int i=hist.size()-1; i>=0; i--)
    Rcpp::Rcout << i+1 << "\t" << hist[i] << "\n";
    
  return(0);
}
