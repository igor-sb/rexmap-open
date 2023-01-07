#include "AlignmentGrid.h"
#include "MergedAlignment.h"
#include "QualityMergingMap.h"

PairedSeqAndQualChar locate_sequence_and_quality_char(
    unsigned int row,
    unsigned int column,
    PairedString &sequences,
    PairedString &qualities
) {
  return PairedSeqAndQualChar(
    sequences.forward[column],
    sequences.reverse[row],
    qualities.forward[column],
    qualities.reverse[row]
  );
}

bool is_path_start(unsigned int &row, unsigned int &column) {
  return column == 0 && row == 0;
}

// [[Rcpp::export]]
std::vector<std::vector<int>> create_scoring_matrix(
    Rcpp::IntegerVector &alignment_scores
) {
  int match = alignment_scores["match"];
  int mismatch = alignment_scores["mismatch"];
  std::vector<std::vector<int>> score_matrix(5, std::vector<int>(5));
  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < 5; j++) {
      score_matrix[i][j] = (i == j && i < 4) ? match : mismatch; 
    }
  }
  return score_matrix;
}

MergedAlignment merge_paired_sequence_and_quality(
    PairedString &paired_sequence,
    PairedString &paired_quality,
    std::vector<int> &path,
    QualityMergingMap quality_merging_map
) {
  unsigned int ncol = paired_sequence.forward.size() + 1;
  unsigned int nrow = paired_sequence.reverse.size() + 1;
  unsigned int flat_id, column = ncol - 1, row = nrow - 1;
  PairedSeqAndQualChar unmerged_pair;
  SeqAndQualChar merged_char;
  MergedAlignment aligned;
  
  while (!is_path_start(row, column)) {
    flat_id = flatten_index(row, column, ncol);
    switch(path[flat_id]) {
    case int('d'): {
      row--; column--;
      unmerged_pair = locate_sequence_and_quality_char(
        row, column, paired_sequence, paired_quality
      );
      merged_char = unmerged_pair.merge(quality_merging_map);
      aligned.append(merged_char).increment_overlap_length();
      if (unmerged_pair.is_equal_seq) aligned.increment_overlap_matches();
      break;
    }
    case int('l'): {
      column--;
      aligned
        .append_sequence(paired_sequence.forward[column])
        .append_quality(paired_quality.forward[column]);
      if (row > 0) aligned.increment_overlap_length();
      break;
    }
    case int('u'): {
      row--;
      aligned
        .append_sequence(paired_sequence.reverse[row])
        .append_quality(paired_quality.reverse[row]);
      if (column < ncol - 1) aligned.increment_overlap_length();
      break;
    }
    default:
      Rcpp::stop("Invalid backtracking value.");
    }
  }
  return aligned.reverse();
};

// [[Rcpp::export]]
Rcpp::List c_align_paired_sequence_and_quality(
    Rcpp::List &sequences,
    Rcpp::List &qualities,
    Rcpp::IntegerVector &alignment_scores,
    std::vector<std::vector<unsigned int>> &match_quality_merging_map,
    std::vector<std::vector<unsigned int>> &mismatch_quality_merging_map
) {
  PairedString sequence_pair(sequences["forward"], sequences["reverse"]);
  PairedString quality_pair(qualities["forward"], qualities["reverse"]);

  AlignmentGrid alignment_grid(
      sequence_pair.forward.size() + 1,
      sequence_pair.reverse.size() + 1,
      "overlap"
  );
  std::vector<std::vector<int>> scoring_matrix = create_scoring_matrix(
    alignment_scores
  );
  int gap_penalty = alignment_scores["gap_penalty"];
  alignment_grid.fill_full_grid(sequence_pair, scoring_matrix, gap_penalty);
  std::vector<int> alignment_grid_path = alignment_grid.get_path();
  
  QualityMergingMap quality_merging_map(
    match_quality_merging_map, mismatch_quality_merging_map
  );
  MergedAlignment alignment = merge_paired_sequence_and_quality(
    sequence_pair, quality_pair, alignment_grid_path, quality_merging_map
  );
  return alignment.as_rcpp_list();
}
