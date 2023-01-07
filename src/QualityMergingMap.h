#pragma once
#include <vector>

class QualityMergingMap {
private:
  
  std::vector<std::vector<unsigned int>> &match;
  std::vector<std::vector<unsigned int>> &mismatch;
  int phred_offset;
  
public:
  
  QualityMergingMap(
    std::vector<std::vector<unsigned int>> &match,
    std::vector<std::vector<unsigned int>> &mismatch,
    int phred_offset = 33
  );
  std::vector<std::vector<unsigned int>> &get_map(bool is_match);
  unsigned int char_to_int(char quality_char);
  char merge_quality_chars(
      char &quality_forward_char,
      char &quality_reverse_char,
      bool is_match
  );
};