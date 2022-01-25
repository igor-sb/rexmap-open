#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

// Nucleotide-to-integer conversion map
static std::map<char, int> nt2int = {{'A',0}, {'C', 1}, {'G', 2}, {'T', 3}, {'N', 4}, {'-', 5}};

// [[Rcpp::export]]
int hamming (std::string s1, std::string s2) {
  // Calculate hamming distance between two sequences of equal length
  int s1_i, s2_i;
  int mismatch = 0;
  size_t len = s1.length();
  // Length check
  if (len != s2.length()) {
    Rprintf("Error: aligned strings are not the same length.\n");
    return 0;
  }
  // Main loop with comparisons
  for (size_t i=0; i<len; i++) {
    s1_i = nt2int[s1.at(i)];
    s2_i = nt2int[s2.at(i)];
    if (s1_i != s2_i) mismatch++;
  }
  return mismatch;
}

// [[Rcpp::export]]
IntegerVector compare_alignment (std::string s1, std::string s2) {
  // alignment strings should not have - in both sequences (like from multiple aln)
  // in the same position.
  int s1_i, s2_i, s1_p, s2_p;
  int match = 0;
  int mismatch = 0;
  int gapopen = 0;
  int gapextend = 0;
  size_t len = s1.length();
  // Length check
  if (len != s2.length()) {
    Rprintf("Error: aligned strings are not the same length.\n");
    return R_NilValue;
  }
  // Keep track of previous symbols for gap extension.
  // Initialize them to any other than 5.
  s1_p = 0;
  s2_p = 0;
  // Main loop with comparisons
  for (size_t i=0; i<len; i++) {
    s1_i = nt2int[s1.at(i)];
    s2_i = nt2int[s2.at(i)];
    if ((s1_i == s2_i) && (s1_i != 4)) { // match except N-N
      match++;
    } else if (s1_i == 5) { // letter is -
      if (s1_p == 5) gapextend++;
      else gapopen++;
    } else if (s2_i == 5) {
      if (s2_p == 5) gapextend++;
      else gapopen++;
    } else { // the only thing left is a mismatch
      mismatch++;
    }
    s1_p = s1_i;
    s2_p = s2_i;
  }
  return IntegerVector::create(match, mismatch, gapopen, gapextend);
}

