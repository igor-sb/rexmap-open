#include "MergedAlignment.h"

const std::string MergedAlignment::get_sequence() {
  return sequence;
}

std::string MergedAlignment::get_quality() {
  return quality;
}

unsigned int MergedAlignment::get_overlap_length() {
  return overlap_length;
}

unsigned int MergedAlignment::get_overlap_matches() {
  return overlap_matches;
}

MergedAlignment &MergedAlignment::append(SeqAndQualChar merged_char) {
  sequence.push_back(merged_char.sequence);
  quality.push_back(merged_char.quality);
  return *this;
}

MergedAlignment &MergedAlignment::append_sequence(char sequence_char) {
  sequence.push_back(sequence_char);
  return *this;
}

MergedAlignment &MergedAlignment::append_quality(char quality_char) {
  quality.push_back(quality_char);
  return *this;
}

MergedAlignment &MergedAlignment::increment_overlap_length() {
  overlap_length++;
  return *this;
}

MergedAlignment &MergedAlignment::increment_overlap_matches() {
  overlap_matches++;
  return *this;
}

MergedAlignment &MergedAlignment::reverse() {
  std::reverse(sequence.begin(), sequence.end());
  std::reverse(quality.begin(), quality.end());
  return *this;
}

Rcpp::List MergedAlignment::as_rcpp_list() {
  return Rcpp::List::create(
    Rcpp::Named("sequence") = get_sequence(),
    Rcpp::Named("quality") = get_quality(),
    Rcpp::Named("overlap_length") = get_overlap_length(),
    Rcpp::Named("overlap_matches") = get_overlap_matches()
  );
}