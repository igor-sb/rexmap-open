#include "mergepairs_helper.h"

// [[Rcpp::export]]
Hashmap<std::vector<int>> test_calc_score_path_first_row(
  std::vector<int> &score,
  std::vector<int> &path,
  int ncol
) {
  calc_score_path_first_row(score, path, ncol);
  return {{"score", score}, {"path", path}};
}

// [[Rcpp::export]]
Hashmap<std::vector<int>> test_calc_score_path_first_column(
    std::vector<int> &score,
    std::vector<int> &path,
    int nrow,
    int ncol,
    int gap_p
) {
    calc_score_path_first_column(score, path, nrow, ncol, gap_p);
    return {{"score", score}, {"path", path}};
}

// [[Rcpp::export]]
Hashmap<std::vector<int>> test_calc_score_path_other(
    std::vector<int> &score,
    std::vector<int> &path,
    Rcpp::List &sequences,
    Vector2d<int> &align_scores,
    int gap_p
) {
  PairedString sequences_pair(sequences["forward"], sequences["reverse"]);
  calc_score_path_other(score, path, sequences_pair, align_scores, gap_p);
  return {{"score", score}, {"path", path}};  
}

// [[Rcpp::export]]
Hashmap<std::vector<int>> test_find_best_scoring_path(
    Rcpp::List &sequences,
    Vector2d<int> scoring_matrix,
    int gap_p
) {
  PairedString sequences_pair(sequences["forward"], sequences["reverse"]);
  return find_best_scoring_path(sequences_pair, scoring_matrix, gap_p);  
}

// [[Rcpp::export]]
MergedAlignment test_merge_by_path_backtrack(
    std::vector<int> &path,
    Rcpp::List &sequences,
    Rcpp::List &qualities,
    Vector2d<unsigned int> &merged_qualities_match,
    Vector2d<unsigned int> &merged_qualities_mismatch
) {
  PairedString sequences_pair(sequences["forward"], sequences["reverse"]);
  PairedString qualities_pair(qualities["forward"], qualities["reverse"]);
  return merge_by_path_backtrack(
    path,
    sequences_pair,
    qualities_pair,
    merged_qualities_match,
    merged_qualities_mismatch
  );
}
