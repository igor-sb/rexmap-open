#include "rexmap.h"
using namespace std;

// [[Rcpp::export]]
void string_compare(string &s1, string &s2) {
  unsigned int count_equal = 0;
  unsigned int count_diff = 0;
  for (size_t i = 0; i < s1.size(); i++) {
    if (s1[i] == s2[i]) {
      count_equal++;
    } else {
      count_diff++;
    }
  }
  return;
}

// [[Rcpp::export]]
void integer_vec_compare(vector<int> &v1, vector<int> &v2) {
  unsigned int count_equal = 0;
  unsigned int count_diff = 0;
  for (size_t i = 0; i < v1.size(); i++) {
    if (v1[i] == v2[i]) {
      count_equal++;
    } else {
      count_diff++;
    }
  }
  return;
}

// [[Rcpp::export]]
void string_vec_compare(vector<string> &v1, vector<string> &v2) {
  unsigned int count_equal = 0;
  unsigned int count_diff = 0;
  for (size_t i = 0; i < v1.size(); i++) {
    if (v1[i] == v2[i]) {
      count_equal++;
    } else {
      count_diff++;
    }
  }
  return;
}

// // [[Rcpp::export]]
// unordered_map<char, unordered_map<char, int>> create_scoring_matrix(
//     int match, int mismatch
// ) {
//   vector<char> letters = {'A', 'C', 'G', 'T', '-'};
//   unordered_map<char, unordered_map<char, int>> scores;
//   for (unsigned int i = 0; i < 5; i++) {
//     for (unsigned int j = 0; j < 5; j++) {
//       scores[letters[i]][letters[j]] = (i == j && i < 4) ? match : mismatch;
//     }
//   }
//   return scores;
// }
