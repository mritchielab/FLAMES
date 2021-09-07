/*
  this file is just for testing, don't even worry about it
  delete it when the time comes
*/
#include "bam.h"
#include <iostream>
#include <htslib/bgzf.h>
#include <htslib/sam.h>
#include <htslib/hts.h>

int main ()
{
  bamFile bam = bam_open("test.bam", "r");

  std::cout << "hello\n";
  return 0;
}