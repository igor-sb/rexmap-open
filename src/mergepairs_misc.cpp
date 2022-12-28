#include "mergepairs.h"
#include <cstdio>
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
int to_int(char &nt) {
  int int_nt;
  switch(nt) {
  case 'A':
    int_nt = 0;
    break;
  case 'C':
    int_nt = 1;
    break;
  case 'G':
    int_nt = 2;
    break;
  case 'T':
    int_nt = 3;
    break;
  default:
    int_nt = 4;
  }
  return(int_nt);
}

// [[Rcpp::export]]
std::vector< std::vector<int> > create_scoring_matrix(int match, int mismatch) {
  std::vector< std::vector<int> > score_matrix(5, std::vector<int>(5));
  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < 5; j++) {
      score_matrix[i][j] = (i == j && i < 4) ? match : mismatch; 
    }
  }
  return score_matrix;
}

// From https://stackoverflow.com/a/18703743/1272087
// [[Rcpp::export]]
std::string join_to_string(std::vector<std::string> str_vector) {
  std::string str;
  for (std::vector<std::string>::const_iterator i = str_vector.begin();
       i != str_vector.end();
       ++i) {
    str += *i;
  }
  return str;  
}

// template <typename T>
void print_flat_vector(const std::vector<int> &v, unsigned int ncol) {
  unsigned int flat_id;
  unsigned int flat_size = v.size();
  char buffer[5];
  
  for (flat_id = 0; flat_id < v.size(); flat_id++) {
    if (flat_id % ncol == 0) std::cout << "\n";
    sprintf(buffer, "% 3d", v[flat_id]);
    std::cout << buffer << " ";
  }
  std::cout << "\n";
  return;
}

// Sum matching integer values from two hashmaps
std::unordered_map<std::string, int> sum_values(
    std::unordered_map<std::string, int> map1,
    std::unordered_map<std::string, int> map2
) {
  std::unordered_map<std::string, int> sum;
  
  for (const auto &key_value : map1) {
    sum[key_value.first] = key_value.second;
  }
  
  for (const auto &key_value : map2) {
    sum[key_value.first] += key_value.second;
  }
  
  return sum;
}
