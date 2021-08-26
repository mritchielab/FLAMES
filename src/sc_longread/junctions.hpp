/*
  Anything related to junctions are goes in here
*/

#include <stdio.h>
#include <iostream>
#include <map>
#include <vector>
#include <any>

typedef struct {
  /*
    a struct used for holding the junctions information
  */
  int left;
  std::vector<int> junctions;
  int right;
} Junctions;

int 
take_closest (std::vector<int> list, int num)
{
  /*
    returns the value in list that is closest to num

    optimize this further later - currently just O(n)
  */
 int output = list.back();
 list.pop_back();

 for (int i : list) 
 {
   if (abs(i - num) < abs(output - num))
   {
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

  output.left = blocks.front().first;
  output.right = blocks.back().second;

  if (blocks.size() > 1)
  {
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
  std::map<std::string, std::vector<int>> all_site;

  for (std::string t : tr_list) {
    if (all_site["left"].size() > 0)
    {
      if (abs(take_closest(all_site["left"], transcript_to_junctions[t].left) - transcript_to_junctions[t].left) > 5)
      {
        all_site["left"].push_back(transcript_to_junctions[t].left);
      }
    }
    else
    {
      all_site["left"].push_back(transcript_to_junctions[t].left);
    }

    if (all_site["right"].size() > 0)
    {
      if (abs(take_closest(all_site["right"], transcript_to_junctions[t].right) - transcript_to_junctions[t].right) > 5)
      {
        all_site["right"].push_back(transcript_to_junctions[t].right);
      }
    }
    else
    {
      all_site["right"].push_back(transcript_to_junctions[t].right);
    }
  }

  Junctions junc;
  return junc;
}