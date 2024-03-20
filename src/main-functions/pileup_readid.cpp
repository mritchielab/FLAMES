#include "htslib/sam.h"
#include <Rcpp.h>
#include <array>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <unistd.h>
#include <unordered_map>
#include <vector>
// [[Rcpp::plugins(cpp17)]]

const std::vector<std::string> BASES = {"A", "T", "C", "G", "del"};

typedef struct plpconf {
  const char *inname;
  samFile *infile;
  sam_hdr_t *in_samhdr;
  hts_idx_t *in_idx;
  hts_itr_t *iter;
} plpconf;

// copied from htslib/samples/pileup.c
// not sure what these are for
inline int plpconstructor(void *data, const bam1_t *b, bam_pileup_cd *cd) {
  return 0;
}
inline int plpdestructor(void *data, const bam1_t *b, bam_pileup_cd *cd) {
  return 0;
}

int readdata(void *data, bam1_t *b) {
  plpconf *conf = (plpconf *)data;
  if (!conf || !conf->infile) {
    return -2; // cant read data
  }

  bool skip = true;
  int ret = 0;
  do {
    ret = conf->iter ? sam_itr_next(conf->infile, conf->iter, b)
                     : sam_read1(conf->infile, conf->in_samhdr, b);
    if (b->core.tid < 0 ||
        (b->core.flag & BAM_FUNMAP)) { // exclude unmapped reads
      skip = true;
      continue;
    } else {
      skip = false;
    }
  } while (skip);

  return ret;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix variant_count_matrix_cpp(Rcpp::String bam_path,
                                             Rcpp::String seqname, int pos,
                                             bool indel,
                                             Rcpp::StringVector barcodes,
                                             bool verbose) {
  if (pos > 1) {
    // htslib is 0-based
    pos = pos - 1;
  }
  plpconf conf = {0};
  bam_plp_t plpiter = NULL;
  const bam_pileup1_t *plp = NULL;

  conf.inname = bam_path.get_cstring();
  const char *seqname_c = seqname.get_cstring();

  if (!(conf.infile = sam_open(conf.inname, "r"))) {
    Rcpp::stop("Failed to open file %s\n", conf.inname);
  }
  if (!(conf.in_samhdr = sam_hdr_read(conf.infile))) {
    Rcpp::stop("Failed to read header from file!\n");
  }
  if (!(plpiter = bam_plp_init(readdata, &conf))) {
    Rcpp::stop("Failed to initialize pileup data\n");
  }
  bam_plp_set_maxcnt(plpiter, INT_MAX); // caps at 2b
  if (!(conf.in_idx = sam_index_load(conf.infile, conf.inname))) {
    Rcpp::stop("Failed to load index for %s", conf.inname);
  }
  int tid = bam_name2id(conf.in_samhdr, seqname_c);
  if ((conf.iter = sam_itr_queryi(conf.in_idx, tid, pos, pos + 1)) == 0) {
    Rcpp::stop("Failed to parse region");
  }

  // set constructor destructor callbacks
  bam_plp_constructor(plpiter, plpconstructor);
  bam_plp_destructor(plpiter, plpdestructor);

  // char seq_nt16_str[] = "=ACMGRSVTWYHKDBN";
  const std::string seq_nt16_string[] = {"=", "A", "C", "M", "G", "R",
                                         "S", "V", "T", "W", "Y", "H",
                                         "K", "D", "B", "N"};

  // fetch an alignment to determine the UMI length first
  bam1_t *bam_tmp = NULL;
  if (!(bam_tmp = bam_init1())) {
    Rcpp::stop("Failed to initialize bamdata");
  }
  if (sam_itr_next(conf.infile, sam_itr_queryi(conf.in_idx, tid, pos, pos + 1),
                   bam_tmp) < 0) {
    Rcpp::stop("Failed to fetch an alignment");
  }
  std::string read_id = bam_get_qname(bam_tmp);
  const std::size_t umi_idx = read_id.find("#");
  const std::size_t id_idx = read_id.find("_");
  if (umi_idx == std::string::npos || id_idx == std::string::npos) {
    Rcpp::stop("Unexpected read id format: %s", read_id);
  }
  if (verbose) {
    Rcpp::Rcout << "Checking read ID format:\n"
                << "ReadID: " << read_id << "\n"
                << "barcode: " << read_id.substr(0, id_idx) << "\n"
                << "UMI: " << read_id.substr(id_idx + 1, umi_idx - id_idx - 1)
                << "\n";
  }
  bam_destroy1(bam_tmp);

  std::unordered_map<std::string, std::unordered_map<std::string, unsigned int>>
      variants;
  int actual_tid = -1, n = -1, actual_pos = -1;
  while ((plp = bam_plp_auto(plpiter, &actual_tid, &actual_pos, &n))) {
    if (actual_tid != tid || actual_pos != pos) {
      continue;
    }

    // iterate all reads and print the read id and the allele
    for (int j = 0; j < n; j++) {

      if (plp[j].is_refskip) {
        continue;
      }

      const std::string read_id = bam_get_qname(plp[j].b);
      const std::string barcode = read_id.substr(0, id_idx);
      const std::string umi = read_id.substr(id_idx + 1, umi_idx - id_idx - 1);

      if (!indel) { // count SNV
        if (plp[j].is_del) {
          variants["del"][barcode]++;
        } else {
          variants[seq_nt16_string[bam_seqi(bam_get_seq(plp[j].b),
                                            plp[j].qpos)]][barcode]++;
        }
      } else { // count indel
        if (plp[j].indel == 0) {
          variants["."][barcode]++;
        } else if (plp[j].indel > 0) {
          // key as '+[inserted_bases]' e.g. '+A' or '+AT'
          std::string inserted_bases = "";
          for (int k = 1; k <= plp[j].indel; k++) {
            inserted_bases += seq_nt16_string[bam_seqi(bam_get_seq(plp[j].b),
                                                       plp[j].qpos + k)];
          }
          variants["+" + inserted_bases][barcode]++;
        } else if (plp[j].indel < 0) {
          // key as '-[deleted_number_of_bases]' e.g. '-2'
          variants["-" + std::to_string(-plp[j].indel)][barcode]++;
        }
      }
    }
  }

  std::vector<std::string> rownames;
  // for SNV, produce consistent rowname order
  // for indel, first row is always "."
  if (!indel) {
    rownames = BASES;
  } else {
    rownames = {"."};
    for (const auto &variant : variants) {
      if (std::find(rownames.begin(), rownames.end(), variant.first) ==
          rownames.end()) {
        rownames.push_back(variant.first);
      }
    }
  }
  Rcpp::NumericMatrix ret(rownames.size(), barcodes.size());
  Rcpp::rownames(ret) = Rcpp::wrap(rownames);
  Rcpp::colnames(ret) = Rcpp::wrap(barcodes);
  for (unsigned int i = 0; i < rownames.size(); i++) {
    for (unsigned int j = 0; j < barcodes.size(); j++) {
      // avoid inserting 0s with variants[rownames[i]][barcodes[j]]
      auto itr = variants[rownames[i]].find(
          (std::string)((Rcpp::String)barcodes[j]).get_cstring());
      if (itr != variants[rownames[i]].end()) {
        ret(i, j) = itr->second;
      }
      // default to 0, no need to update
    }
  }

  if (conf.in_samhdr) {
    sam_hdr_destroy(conf.in_samhdr);
  }
  if (conf.infile) {
    sam_close(conf.infile);
  }
  if (conf.in_idx) {
    hts_idx_destroy(conf.in_idx);
  }
  if (conf.iter) {
    hts_itr_destroy(conf.iter);
  }
  if (plpiter) {
    bam_plp_destroy(plpiter);
  }

  return ret;
}
