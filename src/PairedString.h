#pragma once
#include <string>
/*
 * Container class for forward and reverse DNA sequences and quality scores
 */
class PairedString {
public:
  std::string forward = "";
  std::string reverse = "";
  PairedString() {}
  PairedString(std::string fwd, std::string rev) : forward(fwd), reverse(rev) {}
};