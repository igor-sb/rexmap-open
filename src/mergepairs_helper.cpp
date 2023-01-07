#include "mergepairs_helper.h"
#include "AlignmentGrid.h"

template<typename T>
using Hashmap = std::unordered_map<std::string, T>;

// [[Rcpp::export]]
Hashmap<unsigned int> c_create_flat_index(
  unsigned int row,
  unsigned int column,
  unsigned int num_columns
) {
  Index flat_id(row, column, num_columns);
  return {
    {"current", flat_id.current},
    {"left", flat_id.left},
    {"upper", flat_id.upper},
    {"diag", flat_id.diag}
  };
}

// [[Rcpp::export]]
Hashmap<std::vector<int>> c_alignmentgrid_fill_first_row(
  std::vector<int> &score,
  std::vector<int> &path,
  unsigned int num_rows,
  unsigned int num_columns
) {
  AlignmentGrid score_and_path(num_rows, num_columns, "overlap");
  score_and_path
    .set_manual_score(score)
    .set_manual_path(path)
    .fill_first_row();
  return {
    {"score", score_and_path.get_score()},
    {"path", score_and_path.get_path()}
  };
}

// [[Rcpp::export]]
Hashmap<std::vector<int>> c_alignmentgrid_fill_first_column(
    std::vector<int> &score,
    std::vector<int> &path,
    int num_rows,
    int num_columns,
    int gap_penalty
) {
  AlignmentGrid score_and_path(num_rows, num_columns, "overlap");
  score_and_path
    .set_manual_score(score)
    .set_manual_path(path)
    .fill_first_column(gap_penalty);
  return {
    {"score", score_and_path.get_score()},
    {"path", score_and_path.get_path()}
  };
}

// [[Rcpp::export]]
Hashmap<std::vector<int>> c_alignmentgrid_fill_other(
    std::vector<int> &score,
    std::vector<int> &path,
    Rcpp::List &sequences,
    std::vector<std::vector<int>> &scoring_matrix,
    int gap_penalty
) {
  std::string sequence_fwd = sequences["forward"];
  std::string sequence_rev = sequences["reverse"];
  unsigned int num_columns = sequence_fwd.size() + 1;
  unsigned int num_rows = sequence_rev.size() + 1;
  PairedString sequences_pair(sequence_fwd, sequence_rev);
  AlignmentGrid score_and_path(num_rows, num_columns, "overlap");
  score_and_path
    .set_manual_score(score)
    .set_manual_path(path)
    .fill_other_rows_cols(sequences_pair, scoring_matrix, gap_penalty);
  
  return {
    {"score", score_and_path.get_score()},
    {"path", score_and_path.get_path()}
  };
}

// [[Rcpp::export]]
Hashmap<std::vector<int>> c_alignmentgrid_best_path(
    Rcpp::List &sequences,
    std::vector<std::vector<int>> &scoring_matrix,
    int gap_penalty
) {
  std::string sequence_fwd = sequences["forward"];
  std::string sequence_rev = sequences["reverse"];
  unsigned int num_columns = sequence_fwd.size() + 1;
  unsigned int num_rows = sequence_rev.size() + 1;
  PairedString sequences_pair(sequences["forward"], sequences["reverse"]);
  AlignmentGrid score_and_path(num_rows, num_columns, "overlap");
  score_and_path.find_best_scoring_path(
    sequences_pair, scoring_matrix, gap_penalty
  );
  return {
    {"score", score_and_path.get_score()},
    {"path", score_and_path.get_path()}
  };
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
