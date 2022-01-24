#include "match_cell_barcode.h"

#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <vector>
#include <utility>
#include <unordered_map>
#include <algorithm>
#include <string>
#include <cassert>
#include <fstream>

#include "../utility/edit_dist.h"
#include "../utility/ssw/ssw_cpp.h"
#include "../utility/fastq_utils.h"

using namespace Rcpp;

//static const int MAX_DIST = 2;
static const int BC_LEN = 16;
//static const int UMI_LEN = 10;
static const int BC_WINDOW = 10;

std::string join_path(const std::string p1, const std::string p2)
{
  auto sep = '/';
  return p1.back() == sep ? p1 + p2 : p1 + sep + p2;
}

int find_polyT(std::string &seq, int start_pos)
{
  int i = start_pos;
  while (i < (int)seq.size() - 10)
  {
    if (seq[i] == 'T')
    {
      i++;
    }
    else if (seq.substr(i + 1, 3) == "TTT")
    {
      i = i + 4;
    }
    else if (seq.substr(i + 2, 5) == "TTTTT")
    {
      i = i + 7;
    }
    else
    {
      if (i < (int)seq.size())
      {
        if (i < (start_pos + 40))
        {
          return i;
        }
        else
        {
          if (start_pos + 40 < (int)seq.size())
          {
            return start_pos + 40;
          }
          else
          {
            return start_pos;
          }
        }
      }
      else
      {
        return start_pos;
      }
    }
  }
  return start_pos;
}

char complement(char n)
{
  switch (n)
  {
  case 'A':
    return 'T';
  case 'T':
    return 'A';
  case 'G':
    return 'C';
  case 'C':
    return 'G';
  case 'N':
    return 'N';
  }
  assert(false);
  return ' ';
}

void rc(std::string &seq)
{
  transform(
      begin(seq),
      end(seq),
      begin(seq),
      complement);
  reverse(seq.begin(), seq.end());
}

//int getdir (const char dir[], std::vector<std::string> &files)
//{
//    DIR *dp;
//    struct dirent *dirp;
//    std::string f_name;
//    if((dp  = opendir(dir)) == NULL) {
//        std::cout << "Error(" << errno << ") opening " << dir << std::endl;
//        return errno;
//    }
//
//    while ((dirp = readdir(dp)) != NULL) {
//        f_name = std::string(dirp->d_name);
//        if (f_name[0] != '.') // ingore system files start with `.`
//        {
//          files.push_back(f_name);
//        }
//    }
//    closedir(dp);
//    return 0;
//}

int getdir(const char dir[], std::vector<std::string> &files)
{
  DIR *dp;
  std::string fq_ext1 = "fastq";
  std::string fq_ext2 = "fq";
  struct dirent *dirp;
  std::string f_name;
  if ((dp = opendir(dir)) == NULL)
  {
    Rprintf("Error(%d) opening %s\n", errno, dir);
    //std::cout << "Error(" << errno << ") opening " << dir << std::endl;
    return errno;
  }

  while ((dirp = readdir(dp)) != NULL)
  {
    f_name = std::string(dirp->d_name);
    if (f_name[0] != '.') // ingore system files start with `.`
    {
      std::size_t found1 = f_name.find(fq_ext1);
      std::size_t found2 = f_name.find(fq_ext2);
      if (found1 != std::string::npos)
      {
        files.push_back(f_name);
      }
      else if (found2 != std::string::npos)
      {
        files.push_back(f_name);
      }
    }
  }
  closedir(dp);
  return 0;
}

void get_bc_anno(std::string fn, std::unordered_map<std::string, std::string> &barcode_dict, std::vector<std::string> &barcode_list)
{
  std::ifstream infile(fn);
  std::string line;
  std::string cell_id;
  std::string barcode;
  bool skip_header = false;
  while (std::getline(infile, line))
  {
    std::stringstream linestream(line);
    if (line.size() < 3)
    {
      continue;
    }
    if (line.find(",") != std::string::npos)
    {

      std::getline(linestream, cell_id, ',');
      std::getline(linestream, barcode, ',');
      if (!skip_header)
      {
        skip_header = true; //skip header
        continue;
      }
      if (barcode.length() > 1)
      {
        barcode_dict[barcode] = cell_id;
        barcode_list.push_back(barcode);
      }
    }
    else
    {
      if (line.find("-") != std::string::npos) // 10X barcode file
      {
        std::getline(linestream, barcode, '-');
        barcode_dict[barcode] = barcode;
        barcode_list.push_back(barcode);
      }
      else
      {
        barcode_dict[line] = line;
        barcode_list.push_back(line);
      }
    }
  }
  if (barcode_list.size() < 5)
  {
    Rprintf("Number of cell barcode smaller than 5.\n");
    //std::cout << "Number of cell barcode smaller than 5."<< std::endl;
  }
  else
  {
    Rprintf("First 5 cell barcode:\n");
    //std::cout << "First 5 cell barcode:"<< std::endl;
    for (int i = 0; i < 5; i++)
    {
      Rprintf("\t%s\n", barcode_list[i].c_str());
      //std::cout << "\t" << barcode_list[i]<< std::endl;
    }
  }
}

// Driver function to sort the vector elements
// by second element of pairs
bool sortbysec_dec(const std::pair<int, int> &a,
                   const std::pair<int, int> &b)
{
  return (a.second > b.second);
}

// insensitive to running time.
std::pair<int, int> get_bc_range(std::string fqn, int max_reads, const std::string lf_seq, int MAX_DIST)
{
  std::string query;
  // Declares a default Aligner
  StripedSmithWaterman::Aligner aligner;
  // Declares a default filter
  StripedSmithWaterman::Filter filter;
  // Declares an alignment that stores the result
  StripedSmithWaterman::Alignment alignment;
  std::unordered_map<int, int> bc_pos_dict;
  std::unordered_map<int, int> bc_pos_rev_dict;
  double total_cnt = 0;
  double found_cnt = 0;
  double found_cnt_rev = 0;
  Rprintf("%s\n", fqn.c_str());
  //std::cout << fqn << std::endl;
  gzFile fn = gzopen(fqn.c_str(), "r");
  kseq_t *seq1;
  seq1 = kseq_init(fn);
  while (kseq_read(seq1) >= 0)
  {
    total_cnt++;
    if (total_cnt > max_reads)
    {
      break;
    }
    aligner.Align(seq1->seq.s, lf_seq.c_str(), lf_seq.size(), filter, &alignment, 15);
    if (alignment.mismatches < MAX_DIST + 1)
    {
      found_cnt++;
      if (bc_pos_dict.find(alignment.query_end) != bc_pos_dict.end())
      {
        bc_pos_dict[alignment.query_end]++;
      }
      else
      {
        bc_pos_dict[alignment.query_end] = 1;
      }
    }
    else
    {
      query = std::string(seq1->seq.s);
      rc(query);
      aligner.Align(query.c_str(), lf_seq.c_str(), lf_seq.size(), filter, &alignment, 15);
      if (alignment.mismatches < MAX_DIST + 1)
      {
        found_cnt_rev++;
        if (bc_pos_rev_dict.find(alignment.query_end) != bc_pos_rev_dict.end())
        {
          bc_pos_rev_dict[alignment.query_end]++;
        }
        else
        {
          bc_pos_rev_dict[alignment.query_end] = 1;
        }
      }
    }
  }
  kseq_destroy(seq1);
  gzclose(fn);

  std::vector<std::pair<int, int>> vect;
  std::vector<std::pair<int, int>> vect_rev;
  for (auto &it : bc_pos_dict)
  {
    vect.push_back(std::make_pair(it.first, it.second));
  }
  for (auto &it : bc_pos_rev_dict)
  {
    vect_rev.push_back(std::make_pair(it.first, it.second));
  }
  sort(vect.begin(), vect.end(), sortbysec_dec);
  sort(vect_rev.begin(), vect_rev.end(), sortbysec_dec);

  for (int ix = 0; ix < 20; ix++)
  {
    Rprintf("forward flanking end: %d\t%d\n", vect[ix].first, vect[ix].second);
    //std::cout << "forward flanking end: " << vect[ix].first << "\t"<<vect[ix].second << std::endl;
  }
  for (int ix = 0; ix < 20; ix++)
  {
    Rprintf("reverse comp flanking end: %d\t%d\n", vect_rev[ix].first, vect_rev[ix].second);
    //std::cout << "reverse comp flanking end: " << vect_rev[ix].first << "\t"<<vect_rev[ix].second << std::endl;
  }

  int peak_bc_start = -1;
  int peak_bc_rev_start = -1;

  if (found_cnt / total_cnt > 0.05)
  {
    peak_bc_start = int(vect[0].first + vect[1].first + vect[2].first) / 3;
  }
  if (found_cnt_rev / total_cnt > 0.05)
  {
    peak_bc_rev_start = int(vect_rev[0].first + vect_rev[1].first + vect_rev[2].first) / 3;
  }

  Rprintf("###total read: %f\n", total_cnt);
  Rprintf("###found flanking region: %f\n", found_cnt);
  Rprintf("###found flanking region(rev): %f\n", found_cnt_rev);
  //std::cout << "###total read: " << total_cnt <<std::endl;
  //std::cout << "###found flanking region: " << found_cnt <<std::endl;
  //std::cout << "###found flanking region(rev): " << found_cnt_rev <<std::endl;

  std::pair<int, int> res;
  res = std::make_pair(peak_bc_start, peak_bc_rev_start);
  return res;
}

int get_clo_idx(const char *seq_ptr, int64_t *al, std::vector<int64_t *> &bc_list_ptr, int max_dist)
{
  int i;
  int dis;
  int idx = -1;
  int sml1_dis = 999;
  int sml2_dis = 999;
  for (i = 0; i < 16; ++i)
  {
    al[i] = std::hash<char>{}(seq_ptr[i]);
  }
  for (i = 0; i < (int)bc_list_ptr.size(); i++)
  {
    dis = scutil::edit_distance1(al, 16, bc_list_ptr[i], 16);
    if (dis <= max_dist)
    {
      if (dis < sml1_dis)
      {
        sml1_dis = dis;
        idx = i;
      }
      else if (sml1_dis <= dis && dis < sml2_dis)
      {
        sml2_dis = dis;
      }
    }
  }

  if (sml1_dis < sml2_dis)
  {
    return idx;
  }
  else
  {
    return -1;
  }
}

int get_hm_idx(std::string &q_seq, std::vector<std::string> &barcode_list, int max_dist)
{
  int i;
  int dis;
  int idx = -1;
  int sml1_dis = 999;
  int sml2_dis = 999;
  for (i = 0; i < (int)barcode_list.size(); i++)
  {
    dis = scutil::hamming_distance(q_seq, barcode_list[i]);
    if (dis <= max_dist)
    {
      if (dis < sml1_dis)
      {
        sml1_dis = dis;
        idx = i;
      }
      else if (sml1_dis <= dis && dis < sml2_dis)
      {
        sml2_dis = dis;
      }
    }
  }

  if (sml1_dis < sml2_dis)
  {
    return idx;
  }
  else
  {
    return -1;
  }
}


// [[Rcpp::export]]
void match_cell_barcode(String fastq_dir, String stats_file, String out_fastq, String ref_csv, int MAX_DIST, int UMI_LEN = 10)
{
  // "usage: <1.fastq folder> <2.output cell barcode statistics file> <3.fastq output reads that matched cell barcode> <4.barcode reference from short read 10X data> <5.max edit distance> [6. UMI length (default: 10)]"

  gzFile o_stream_gz = gzopen(out_fastq.get_cstring(), "wb2"); // open the output file, writing it in binary mode with compression level 2
  int total_cnt = 0;
  int found_cnt = 0;
  int do_match = 0;
  int not_match = 0;
  int do_match_hm = 0;
  int too_short = 0;
  int bc_idx;
  std::string query;
  // Declares a default Aligner
  StripedSmithWaterman::Aligner aligner;
  // Declares a default filter
  StripedSmithWaterman::Filter filter;
  // Declares an alignment that stores the result
  StripedSmithWaterman::Alignment alignment;
  std::unordered_map<std::string, std::string> barcode_dict;
  std::vector<std::string> barcode_list;
  std::string bc_anno_fn = ref_csv.get_cstring();
  get_bc_anno(bc_anno_fn, barcode_dict, barcode_list);

  std::unordered_map<std::string, int> barcode_res;
  std::unordered_map<std::string, std::string> fuzzy_barcode_match;

  std::unordered_map<int, int> polyT_pos_stat;

  const std::string ref_left = "CTACACGACGCTCTTCCGATCT";

  int64_t *al = (int64_t *)malloc(16 * sizeof(int64_t));
  std::vector<int64_t *> bc_list_ptr;

  for (auto &it : barcode_list)
  {
    barcode_res[it] = 0;
    bc_list_ptr.push_back((int64_t *)malloc(16 * sizeof(int64_t)));
    for (int m = 0; m < 16; m++)
    {
      bc_list_ptr.back()[m] = std::hash<char>{}(it[m]);
    }
  }
  std::string bc_string;
  std::string true_bc;
  int matched_idx, polyT_idx;

  std::string name, seq, qual;

  std::vector<std::string> files;
  getdir(fastq_dir.get_cstring(), files);
  bool found_match = false;
  std::string seed_file = files[0];
  for (auto &it : files)
  {
    if (it.find("pass") != std::string::npos)
    {
      seed_file = it;
    }
  }

  std::pair<int, int> bc_range = get_bc_range(join_path(std::string(fastq_dir.get_cstring()), seed_file).c_str(), 30000, ref_left, MAX_DIST);
  //std::cout << bc_range.first <<" @@@@@@ " <<  bc_range.second << std::endl;
  Rprintf("%d @@@@@@ %d\n", bc_range.first, bc_range.second);
  for (auto &it : files)
  {
    Rprintf("%s\n", it.c_str());
    //std::cout << it << std::endl;
    //get_bc_range(join_path(std::string(fastq_dir.get_cstring(),it), 100000, ref_left);
    //break;
    gzFile fn = gzopen(join_path(std::string(fastq_dir.get_cstring()), it).c_str(), "r");
    kseq_t *seq1;
    seq1 = kseq_init(fn);
    while (kseq_read(seq1) >= 0)
    {

      total_cnt++;
      if (total_cnt % 500000 == 0)
      {
        Rprintf("%d :: %d :: %d\n", total_cnt, do_match_hm, do_match);
        //std::cout << total_cnt << " :: " << do_match_hm << " :: " << do_match << std::endl;
      }
      found_match = false;
      seq = std::string(seq1->seq.s);
      qual = std::string(seq1->qual.s);
      name = std::string(seq1->name.s);

      if ((int)seq.size() < std::max(bc_range.first, bc_range.second) + BC_LEN + UMI_LEN + 20)
      {
        too_short++;
        continue;
      }

      if (bc_range.first > -1)
      {
        for (bc_idx = std::max(bc_range.first - BC_WINDOW, 0); bc_idx < std::min(bc_range.first + BC_WINDOW, (int)seq.size() - BC_LEN); bc_idx++)
        {
          bc_string = seq.substr(bc_idx, BC_LEN);
          if (fuzzy_barcode_match.find(bc_string) != fuzzy_barcode_match.end())
          {
            barcode_res[fuzzy_barcode_match[bc_string]]++;
            true_bc = fuzzy_barcode_match[bc_string];
            do_match_hm++;
            found_match = true;
            break;
          }
          else
          {
            matched_idx = get_hm_idx(bc_string, barcode_list, MAX_DIST);
            if (matched_idx > 0)
            {
              true_bc = barcode_list[matched_idx];
              do_match_hm++;
              found_match = true;
              barcode_res[barcode_list[matched_idx]]++;
              fuzzy_barcode_match[bc_string] = barcode_list[matched_idx];
              break;
            }
          }
        }
      }

      if (found_match)
      {
        name = true_bc + "_" + seq.substr(bc_idx + BC_LEN, UMI_LEN) + "#" + name;
        polyT_idx = find_polyT(seq, bc_idx + BC_LEN + UMI_LEN + 2);
        //polyT_idx = bc_idx+BC_LEN+UMI_LEN+2;
        polyT_pos_stat[polyT_idx - (bc_idx + BC_LEN + UMI_LEN + 2)]++;
        seq = seq.substr(polyT_idx, seq.size() - polyT_idx);
        qual = qual.substr(polyT_idx, qual.size() - polyT_idx);
        fq_gz_write(o_stream_gz, name, qual, seq);
        continue;
      }

      if (!found_match)
      {
        aligner.Align(seq1->seq.s, ref_left.c_str(), ref_left.size(), filter, &alignment, 15);
        if ((alignment.mismatches < MAX_DIST + 1) & (alignment.query_end + BC_LEN + UMI_LEN + 50 < (int)seq.size()))
        {
          found_cnt++;
          bc_string = std::string((seq1->seq.s) + alignment.query_end, 16);
          if (fuzzy_barcode_match.find(bc_string) != fuzzy_barcode_match.end())
          {
            barcode_res[fuzzy_barcode_match[bc_string]]++;
            true_bc = fuzzy_barcode_match[bc_string];
            do_match++;
            found_match = true;
          }
          else
          {
            matched_idx = get_clo_idx((seq1->seq.s) + alignment.query_end, al, bc_list_ptr, MAX_DIST + 1);
            if (matched_idx > 0)
            {
              true_bc = barcode_list[matched_idx];
              do_match++;
              found_match = true;
              barcode_res[barcode_list[matched_idx]]++;
              fuzzy_barcode_match[bc_string] = barcode_list[matched_idx];
            }
          }
        }
        else
        {
          rc(seq);
          reverse(qual.begin(), qual.end());

          if (bc_range.second > -1) // look at reverse strand
          {

            for (bc_idx = std::max(bc_range.second - 5, 0); bc_idx < std::min(bc_range.second + 5, (int)seq.size() - BC_LEN); bc_idx++)
            {
              bc_string = seq.substr(bc_idx, BC_LEN);
              if (fuzzy_barcode_match.find(seq.substr(bc_idx, BC_LEN)) != fuzzy_barcode_match.end())
              {
                barcode_res[fuzzy_barcode_match[bc_string]]++;
                true_bc = fuzzy_barcode_match[bc_string];
                do_match_hm++;
                found_match = true;
                break;
              }
              else
              {
                matched_idx = get_hm_idx(bc_string, barcode_list, MAX_DIST);
                if (matched_idx > 0)
                {
                  true_bc = barcode_list[matched_idx];
                  do_match_hm++;
                  found_match = true;
                  barcode_res[barcode_list[matched_idx]]++;
                  fuzzy_barcode_match[bc_string] = barcode_list[matched_idx];
                  break;
                }
              }
            }
            if (found_match)
            {
              name = true_bc + "_" + seq.substr(bc_idx + BC_LEN, UMI_LEN) + "#" + name;
              polyT_idx = find_polyT(seq, bc_idx + BC_LEN + UMI_LEN + 2);
              //polyT_idx = bc_idx+BC_LEN+UMI_LEN+2;
              polyT_pos_stat[polyT_idx - (bc_idx + BC_LEN + UMI_LEN + 2)]++;
              seq = seq.substr(polyT_idx, seq.size() - polyT_idx);
              qual = qual.substr(polyT_idx, qual.size() - polyT_idx);
              fq_gz_write(o_stream_gz, name, qual, seq);
              continue;
            }
          }

          if (!found_match)
          {
            aligner.Align(seq.c_str(), ref_left.c_str(), ref_left.size(), filter, &alignment, 15);
            if ((alignment.mismatches < MAX_DIST + 1) & (alignment.query_end + BC_LEN + UMI_LEN + 50 < (int)seq.size()))
            {
              found_cnt++;
              bc_string = seq.substr(alignment.query_end, 16);
              if (fuzzy_barcode_match.find(bc_string) != fuzzy_barcode_match.end())
              {
                barcode_res[fuzzy_barcode_match[bc_string]]++;
                true_bc = fuzzy_barcode_match[bc_string];
                do_match++;
                found_match = true;
              }
              else
              {
                matched_idx = get_clo_idx(bc_string.c_str(), al, bc_list_ptr, MAX_DIST + 1);
                if (matched_idx > 0)
                {
                  true_bc = barcode_list[matched_idx];
                  do_match++;
                  found_match = true;
                  barcode_res[barcode_list[matched_idx]]++;
                  fuzzy_barcode_match[bc_string] = barcode_list[matched_idx];
                }
              }
            }
          }
        }
        if (found_match)
        {
          name = true_bc + "_" + seq.substr(alignment.query_end + BC_LEN, UMI_LEN) + "#" + name;
          polyT_idx = find_polyT(seq, alignment.query_end + BC_LEN + UMI_LEN + 2);
          //polyT_idx = alignment.query_end+BC_LEN+UMI_LEN+2;
          polyT_pos_stat[polyT_idx - (alignment.query_end + BC_LEN + UMI_LEN + 2)]++;
          seq = seq.substr(polyT_idx, seq.size() - polyT_idx);
          qual = qual.substr(polyT_idx, qual.size() - polyT_idx);
          if (seq.size() > 50)
          {
            fq_gz_write(o_stream_gz, name, qual, seq);
          }
        }
        else
        {
          not_match++;
        }
      }
    }
    kseq_destroy(seq1);
    gzclose(fn);
    //break;
  }
  std::ofstream ofile(stats_file.get_cstring());
  ofile << "cell_id,count\n";
  for (auto &it : barcode_res)
  {
    if (it.second > 0)
    {
      ofile << it.first << "," << it.second << "\n";
    }
  }
  ofile.close();
  gzclose(o_stream_gz);
  Rprintf("###polyT length stat: \n");
  //std::cout << "###polyT length stat: " <<std::endl;
  for (auto &it : polyT_pos_stat)
  {
    Rprintf("\t%d\t%d\n", it.first, it.second);
    //std::cout << "\t" << it.first << "\t" << it.second <<std::endl;
  }

  Rprintf("###total read: %d\n", total_cnt);
  Rprintf("###barcode hm match: %d\n", do_match_hm);
  Rprintf("###barcode fuzzy match: %d\n", do_match);
  Rprintf("###barcode not match: %d\n", not_match);
  Rprintf("###too short: %d\n", too_short);

  // free everything when we're done
  free(al);
  for (const auto & entry : bc_list_ptr) {
    free(entry);
  }

  //std::cout << "###total read: " << total_cnt <<std::endl;
  //std::cout << "###barcode hm match: " << do_match_hm <<std::endl;
  //std::cout << "###barcode fuzzy match: " << do_match <<std::endl;
  //std::cout << "###barcode not match: " << not_match <<std::endl;
  //std::cout << "###too short: " << too_short <<std::endl;
}
