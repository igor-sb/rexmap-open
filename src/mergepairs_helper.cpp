#include "mergepairs_helper.h"

// [[Rcpp::export]]
std::unordered_map<std::string, std::vector<int>>
test_calc_score_path_first_row(
  std::vector<int> &score,
  std::vector<int> &path,
  int ncol
) {
  calc_score_path_first_row(score, path, ncol);
  return std::unordered_map<std::string, std::vector<int>> {
    {"score", score},
    {"path", path}
  };
}

// [[Rcpp::export]]
std::unordered_map<std::string, std::vector<int>>
  test_calc_score_path_first_column(
    std::vector<int> score,
    std::vector<int> path,
    int nrow,
    int ncol,
    int gap_p
  ) {
    calc_score_path_first_column(score, path, nrow, ncol, gap_p);
    return std::unordered_map<std::string, std::vector<int>> {
      {"score", score},
      {"path", path}
    };
  }

// [[Rcpp::export]]
std::unordered_map<std::string, std::vector<int>> test_calc_score_path_other(
    std::vector<int> score,
    std::vector<int> path,
    std::vector<std::string> sequences,
    std::vector<std::vector<int>> align_scores,
    int gap_p
) {
  calc_score_path_other(score, path, sequences, align_scores, gap_p);
  return std::unordered_map<std::string, std::vector<int>> {
    {"score", score},
    {"path", path}
  };  
}

