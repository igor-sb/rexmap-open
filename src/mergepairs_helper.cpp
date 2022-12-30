#include "mergepairs_helper.h"

// [[Rcpp::export]]
std::unordered_map<std::string, std::vector<int>>
test_calc_score_path_first_row(
  std::vector<int> &score,
  std::vector<int> &path,
  int ncol
) {
  calc_score_path_first_row(score, path, ncol);
  return {{"score", score}, {"path", path}};
}

// [[Rcpp::export]]
std::unordered_map<std::string, std::vector<int>>
  test_calc_score_path_first_column(
    std::vector<int> &score,
    std::vector<int> &path,
    int nrow,
    int ncol,
    int gap_p
  ) {
    calc_score_path_first_column(score, path, nrow, ncol, gap_p);
    return {{"score", score}, {"path", path}};
  }

std::unordered_map<std::string, std::string> vector_to_unordered_map(
  std::vector<std::string> &sequences
) {
  return {{"forward", sequences[0]}, {"reverse", sequences[1]}};
}

// [[Rcpp::export]]
std::unordered_map<std::string, std::vector<int>> test_calc_score_path_other(
    std::vector<int> &score,
    std::vector<int> &path,
    std::vector<std::string> &sequences,
    std::vector<std::vector<int>> &align_scores,
    int gap_p
) {
  std::unordered_map<std::string, std::string> sequences_umap = 
    vector_to_unordered_map(sequences);
  calc_score_path_other(score, path, sequences_umap, align_scores, gap_p);
  return {{"score", score}, {"path", path}};  
}

// [[Rcpp::export]]
std::unordered_map<std::string, std::vector<int>> test_find_best_scoring_path(
    std::vector<std::string> &sequences,
    std::vector<std::vector<int>> scoring_matrix,
    int gap_p
) {
  std::unordered_map<std::string, std::string> sequences_umap = 
    vector_to_unordered_map(sequences);
  return find_best_scoring_path(sequences_umap, scoring_matrix, gap_p);  
}

// [[Rcpp::export]]
std::unordered_map<std::string, std::string> test_merge_by_path_backtrack(
    std::vector<int> &path,
    std::vector<std::string> &sequences,
    std::vector<std::string> &qualities,
    std::vector<std::vector<unsigned int>> &merged_qualities_match,
    std::vector<std::vector<unsigned int>> &merged_qualities_mismatch
) {
  std::unordered_map<std::string, std::string> sequences_umap = 
    vector_to_unordered_map(sequences);
  std::unordered_map<std::string, std::string> qualities_umap = 
    vector_to_unordered_map(qualities);
  return merge_by_path_backtrack(
    path,
    sequences_umap,
    qualities_umap,
    merged_qualities_match,
    merged_qualities_mismatch
  );
}
