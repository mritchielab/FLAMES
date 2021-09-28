/*
  Anything related to junctions are goes in here
*/
#include "junctions.h"


int 
take_closest (std::vector<int> list, int num)
{
  /*
    returns the value in list that is closest to num
  */

  int output = list.back();
  list.pop_back();

  for (const auto & i : list) {
    if (abs(i - num) < abs(output - num)) {
      output = i;
    }
  }
  return output;
}


Junctions 
blocks_to_junctions (std::vector<std::pair<int, int>> blocks)
{
  /*
    takes the blocks,
    converts them into a junctions object
  */

  Junctions output;

  output.left = {blocks.front().first};
  output.right = {blocks.back().second};

  if (blocks.size() > 1) {
    for (int i = 1; i < blocks.size(); i++) {
      output.junctions.push_back(blocks[i - 1].second);
      output.junctions.push_back(blocks[i].first);
    }
  }

  return output;
}


Junctions 
get_TSS_TES_site (std::map<std::string, Junctions> transcript_to_junctions, std::vector<std::string> tr_list)
{
  Junctions
  all_site;

  for (const auto & t : tr_list) {
    if (all_site.left.size() > 0) {
      if (abs(take_closest(all_site.left, transcript_to_junctions[t].left[0]) - transcript_to_junctions[t].left[0]) > 5) {
        all_site.left.push_back(transcript_to_junctions[t].left[0]);
      }
    } else {
      all_site.left.push_back(transcript_to_junctions[t].left[0]);
    }

    if (all_site.right.size() > 0) {
      if (abs(take_closest(all_site.right, transcript_to_junctions[t].right[0]) - transcript_to_junctions[t].right[0]) > 5) {
        all_site.right.push_back(transcript_to_junctions[t].right[0]);
      }
    } else {
      all_site.right.push_back(transcript_to_junctions[t].right[0]);
    }
  }

  return all_site;
}


std::set<int> 
get_splice_site (std::map<std::string, Junctions> transcript_to_junctions, std::vector<std::string> tr_list)
{
  std::set<int>
  all_site;

  // add all of the junctions of 
  // everything from tr_list to all_site
  for (std::string t : tr_list) {
    for (const auto & junctions : transcript_to_junctions[t].junctions) {
      all_site.insert(junctions);
    }
  }

  return all_site;
}
