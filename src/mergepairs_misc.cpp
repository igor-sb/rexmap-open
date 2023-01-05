#include <Rcpp.h>
#include "rexmap.h"
#include <cstdio>
// [[Rcpp::plugins(cpp11)]]

// 2D to 1D array index conversion
unsigned int flat_index(
    unsigned int column_index,
    unsigned int row_index,
    unsigned int ncol
) {
  return row_index * ncol + column_index;
}

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
Vector2d<int> create_scoring_matrix(int match, int mismatch) {
  Vector2d<int> score_matrix(5, std::vector<int>(5));
  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < 5; j++) {
      score_matrix[i][j] = (i == j && i < 4) ? match : mismatch; 
    }
  }
  return score_matrix;
}

// template <typename T>
void print_flat_vector(const std::vector<int> &v, unsigned int ncol) {
  unsigned int flat_id;
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
Hashmap<int> sum_values(
    Hashmap<int> map1,
    Hashmap<int> map2
) {
  Hashmap<int> sum;
  
  for (const auto &key_value : map1) {
    sum[key_value.first] = key_value.second;
  }
  
  for (const auto &key_value : map2) {
    sum[key_value.first] += key_value.second;
  }
  
  return sum;
}
