#include "mergepairs_helper.h"
#include "AlignmentGrid.h"
#include <Rcpp.h>

// [[Rcpp::export]]
std::unordered_map<std::string, unsigned int> c_create_flat_index(
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
std::unordered_map<std::string, std::vector<int>>
  c_alignmentgrid_fill_first_row(
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
std::unordered_map<std::string, std::vector<int>>
  c_alignmentgrid_fill_first_column(
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
std::unordered_map<std::string, std::vector<int>> c_alignmentgrid_fill_other(
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
std::unordered_map<std::string, std::vector<int>> c_alignmentgrid_best_path(
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
  score_and_path.fill_full_grid(
    sequences_pair, scoring_matrix, gap_penalty
  );
  return {
    {"score", score_and_path.get_score()},
    {"path", score_and_path.get_path()}
  };
}

// [[Rcpp::export]]
Rcpp::List c_merge_paired_sequence_and_quality(
    Rcpp::List &sequence,
    Rcpp::List &quality,
    std::vector<int> &grid_path,
    std::vector<std::vector<unsigned int>> &match_quality_merging_map,
    std::vector<std::vector<unsigned int>> &mismatch_quality_merging_map
) {
  PairedString sequence_pair(sequence["forward"], sequence["reverse"]);
  PairedString quality_pair(quality["forward"], quality["reverse"]);
  QualityMergingMap quality_merging_map(
      match_quality_merging_map, mismatch_quality_merging_map
  );
  MergedAlignment alignment = merge_paired_sequence_and_quality(
    sequence_pair, quality_pair, grid_path, quality_merging_map
  );
  return alignment.as_rcpp_list();
}
