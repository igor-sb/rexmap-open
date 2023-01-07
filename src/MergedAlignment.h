#pragma once
#include <Rcpp.h>
#include "PairedString.h"
#include "SeqAndQual.h"

unsigned int flatten_index(
    unsigned int row,
    unsigned int column,
    unsigned int num_columns
);

PairedSeqAndQualChar locate_sequence_and_quality_char(
    unsigned int row,
    unsigned int column,
    PairedString &sequences,
    PairedString &qualities
);

class MergedAlignment {
private:
  
  std::string sequence = "";
  std::string quality = "";
  unsigned int overlap_length = 0;
  unsigned int overlap_matches = 0;
  
public:
  
  MergedAlignment() {}
  const std::string get_sequence();
  std::string get_quality();
  unsigned int get_overlap_length();
  unsigned int get_overlap_matches();
  MergedAlignment &append(SeqAndQualChar merged_char);
  MergedAlignment &append_sequence(char sequence_char);
  MergedAlignment &append_quality(char quality_char);
  MergedAlignment &increment_overlap_length();
  MergedAlignment &increment_overlap_matches();
  MergedAlignment &reverse();
  MergedAlignment merge_pairs(
      std::vector<int> &path,
      PairedString &paired_sequence,
      PairedString &paired_quality,
      QualityMergingMap quality_merging_map
  );
  Rcpp::List as_rcpp_list();
};

