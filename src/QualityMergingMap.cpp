#include "QualityMergingMap.h"

QualityMergingMap::QualityMergingMap(
    std::vector<std::vector<unsigned int>> &match,
    std::vector<std::vector<unsigned int>> &mismatch,
    int phred_offset
  ) : match(match), mismatch(mismatch), phred_offset(phred_offset) {}
  
std::vector<std::vector<unsigned int>> &QualityMergingMap::get_map(
    bool is_match
) {
  return is_match ? match : mismatch;
}
  
unsigned int QualityMergingMap::char_to_int(char quality_char) {
  return int(quality_char) - phred_offset;
}
  
char QualityMergingMap::merge_quality_chars(
    char &quality_forward_char,
    char &quality_reverse_char,
    bool is_match
) {
  unsigned int forward_id = char_to_int(quality_forward_char) - 1;
  unsigned int reverse_id = char_to_int(quality_reverse_char) - 1;
  unsigned int larger_id = forward_id > reverse_id ? forward_id : reverse_id;
  unsigned int smaller_id = forward_id > reverse_id ? reverse_id : forward_id;
  return char(get_map(is_match)[larger_id][smaller_id] + phred_offset);
};