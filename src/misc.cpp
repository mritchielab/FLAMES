/*
  miscellaneous helper functions and things
  they will not live here forever, 
  once I figure out how to categorize them better then I will move them
*/
#include "misc.h"

std::vector<int>
find_best_splice_chain(std::vector<int> raw_iso, std::vector<std::vector<int>> junction_list, int MAX_DIST)
{
  int best_match[3] = {-1, 0, 0};
  for (int i = 0; i < junction_list.size(); i++) {
    // hold onto this
    std::vector<int>
    junction = junction_list[i];
    
    // populate vector i_st with the indices of junction_list[i] entries
    // that match certain criteria
    std::vector<int> iter_start;
    for (int j = 0; j < junction_list[i].size(); j++) {

      if (abs((junction[j] - raw_iso[1]) < MAX_DIST)) {
        iter_start.push_back(j);
      }
    }

    if (iter_start.size() == 1) {
      int iter_start_int = iter_start[0];
      int iter = iter_start_int + 1;
      int iter_end = iter_start_int;
      while (iter < junction.size() && iter - iter_start_int + 1 < raw_iso.size()) {
        if (abs(junction[iter] - raw_iso[iter - iter_start_int + 1]) < MAX_DIST) {
          iter_end = iter;
          iter += 1;
        } else {
          break;
        }
      }

      if (iter_end - iter_start_int >= best_match[1]) {
        best_match[0] = i;
        best_match[1] = iter_end - iter_start_int;
        best_match[2] = iter_start_int;
      }
    }
  }

  if (best_match[0] >= 0 && best_match[1] >= 3) {
    std::vector<int>
    updated_iso = raw_iso;

    for (int i = best_match[2]; i < best_match[2] + best_match[1] + 1; i++) {
      updated_iso[i - best_match[2] + 1] = junction_list[best_match[0]][i];
    }

    return updated_iso;
  } else {
    return raw_iso;
  }
}

int
if_exon_contains(std::vector<int> s1, std::vector<int> s2, int MAX_TOLERANCE)
{
  /*
    checks if s2 is in s1
    searching for exact match
  */
  
  if (s2.size() == 2) { // ignore single exon transcripts 
    return 0;
  }

  auto fs_elem = std::find(s1.begin(), s1.end(), s2[1]);

  if (fs_elem == s1.end()) { // ignore if s2[1] is not in s1 
    return 0;
  }

  // get the index of the element
  int fs = std::distance(s1.begin(), fs_elem);

  if ((fs == 0) || ((s2[0] - s1[fs - 1]) < -MAX_TOLERANCE)) { // ignore if left is not within s1 
    return 0;
  }

  for (int i = 2; i < s2.size() - 1; i++) {
    if (fs + i - 1 > s1.size() - 1) {
      return 0;
    }

    if (s1[fs + i - 1] != s2[i]) {
      return 0;
    }
  }

  if (fs + s2.size() - 2 > s1.size() - 1) {
    return 0;
  }
  if ((s2.back() - s1[fs + s2.size() - 2]) > MAX_TOLERANCE) {
    return 0;
  }

  // if all of this is true, return true
  return 1;
}

float
get_exon_sim_pct(std::vector<int> exons1, std::vector<int> exons2)
{
  /*
    takes two exon transcripts,
    returns the percentage of coverage between them
  */

  auto
  pos_overlap = [] (std::pair<int, int> pos1, std::pair<int, int> pos2) 
  {
    if ((pos1.second <= pos2.first) ||
        (pos1.first >= pos2.second)) {
      return 0;
    } else {
      return (std::min(pos1.second, pos2.second) - std::max(pos1.first, pos2.first));
    }
  };

  auto
  sum_of_exon = [] (std::vector<int> exon)
  {
    /* add up the pairs of an exon */

    int sum = 0;
    auto pairs = pairwise(exon);
    for (const auto & pair : pairs) {
      sum += pair.second - pair.first;
    }

    return sum;
  };

  auto e1_len = sum_of_exon(exons1);
  auto e2_len = sum_of_exon(exons2);

  float total = 0;
  for (const auto & pair1 : pairwise(exons1)) {
    for (const auto & pair2 : pairwise(exons2)) {
      total += pos_overlap(pair1, pair2);
    }
  }
  return total / std::max(e1_len, e2_len);
}

std::vector<std::pair<int, int>>
pairwise (std::vector<int> input)
{
  /*
    takes a vector,
    splits it up into a vector of pairs
    {1, 2, 3, 4, 5} -> {{1, 2}, {3, 4}}
  */
  std::vector<std::pair<int, int>>
  output;

  for (int i = 1; i < input.size(); i+=2) {
    std::pair<int, int>
    new_pair = {input[i-1], input[i]};

    output.push_back(new_pair);
  }

  return output;
}

int
iv_overlap (std::pair<int, int> iv1, std::pair<int, int> iv2)
{
  /* takes two ivs as pairs, calculates the overlap between them */

  return std::max(0, std::min(iv1.second, iv2.second) - std::max(iv2.first, iv1.first));
}

int
exon_overlap (std::vector<int> exons1, std::vector<int> exons2)
{
  /* takes two exons, returns the total overlap between them */

  int total = 0;
  for (const auto & e1 : pairwise(exons1)) {
    for (const auto & e2 : pairwise(exons2)) {
      total += iv_overlap(e1, e2);
    } 
  }
  return total;
}

std::vector<std::pair<std::string, std::string>>
get_fa(std::string filename)
{
  auto
  strip = [] (std::string input)
  {
    /* removes leading and trailing whitespace */

    auto start = input.begin();
    auto end = input.rbegin();

    while (isspace(*start)) {++start;}
    while (isspace(*end)) {++end;}

    return std::string(start, end.base());
  };

  auto
  string_toupper = [] (std::string input)
  {
    /* converts a string to uppercase */

    std::string output = input;

    for (int i = 0; i < input.size(); i++) {
      output[i] = std::toupper(output[i]);
    }
    return output;
  };

  std::string ch = "";
  std::vector<std::string> sequence = {};
  std::vector<std::pair<std::string, std::string>> output;

  std::ifstream infile(filename);

  std::string line;
  while (std::getline(infile, line)) {
    if (line[0] == '>') {
      if (ch != "") {
        std::stringstream sequence_string;
        for (const auto & i : sequence) {
          sequence_string << i;
        }

        output.push_back({ch, sequence_string.str()});
      }
      ch = string_toupper(strip(std::string(line.begin()+1, line.end())));
      sequence = {};
    } else {
      sequence.push_back(string_toupper(strip(line)));
    }
  }

  return output;
}