/*
  miscellaneous helper functions and things
*/

#include <string>

// Oliver already wrote this one 
// - it's just here until we merge our branches
struct Pos 
{
  /*
    this is a struct to populate junction_dict
  */
  std::string chr;
  int start, end;
  char strand;
  std::string parent_id;
};


  struct Iso
  {
    /*
      a data container used in Isoforms,
      specifically for known_isoforms and match_known_annotation
    */

    long support_cnt;
    std::string transcript_id; 
    std::string gene_id;

    Iso (long support_cnt_, std::string transcript_id_, std::string gene_id_)
    {
      support_cnt = support_cnt_;
      transcript_id = transcript_id_;
      gene_id = gene_id_;
    };
  };