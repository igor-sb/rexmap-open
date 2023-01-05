#include <Rcpp.h>
#include "rexmap.h"

struct SequenceChar {
  char sequence;
  char quality;
};

struct PairedChar {
  SequenceChar forward;
  SequenceChar reverse;
};

Hashmap<char> get_current_chars(
    unsigned int row,
    unsigned int col,
    PairedString &sequences,
    PairedString &qualities
) {
  return {
    {"sequence_forward", (sequences.forward)[col]},
    {"quality_forward", (qualities.forward)[col]},
    {"sequence_reverse", (sequences.reverse)[row]},
    {"quality_reverse", (qualities.reverse)[row]}
  };
}

// [[Rcpp::export]]
char get_merged_qualities(
    char &q1char,
    char &q2char,
    Vector2d<unsigned int> &merged_qualities
) {
  unsigned int q1 = int(q1char) - PHRED_OFFSET;
  unsigned int q2 = int(q2char) - PHRED_OFFSET;
  unsigned int q_merged;
  if (q2 > q1) std::swap(q1, q2);
  q_merged = merged_qualities[q1 - 1][q2 - 1] + PHRED_OFFSET;
  return char(q_merged);
}

Hashmap<char> merge_forward_and_reverse_chars(
    Hashmap<char> chars,
    Vector2d<unsigned int> merged_qualities_match,
    Vector2d<unsigned int> merged_qualities_mismatch
) {
  Hashmap<char> merged;
  if (chars["sequence_forward"] == chars["sequence_reverse"]) {
    merged["sequence"] = chars["sequence_forward"];
    merged["quality"] = get_merged_qualities(
      chars["quality_forward"],
      chars["quality_reverse"],
      merged_qualities_match
    );
  } else {
    if (chars["quality_forward"] < chars["quality_reverse"]) {
      merged["sequence"] = chars["sequence_reverse"];
    } else {
      merged["sequence"] = chars["sequence_forward"];
    }
    merged["quality"] = get_merged_qualities(
      chars["quality_forward"],
      chars["quality_reverse"],
      merged_qualities_mismatch
    );
  }
  return merged;
}

bool is_first_row_or_last_col(
    unsigned int &col,
    unsigned int &row,
    unsigned int &ncol
) {
  return row == 0 || col == ncol - 1;
}

bool is_same_sequence_char(
    std::unordered_map<std::string, char> &unmerged_chars
) {
  return unmerged_chars["sequence_forward"] == 
    unmerged_chars["sequence_reverse"];
}

bool is_path_start(unsigned int &row, unsigned int &col) {
  return col == 0 && row == 0;
}

MergedAlignment merge_by_path_backtrack(
    std::vector<int> &path,
    PairedString &sequences,
    PairedString &qualities,
    Vector2d<unsigned int> &merged_qualities_match,
    Vector2d<unsigned int> &merged_qualities_mismatch
) {
  unsigned int ncol = sequences.forward.size() + 1;
  unsigned int nrow = sequences.reverse.size() + 1;
  unsigned int flat_id, col = ncol - 1, row = nrow - 1;
  MergedAlignment aligned;
  aligned.overlap_length = 0;
  aligned.overlap_matches = 0;
  Hashmap<char> unmerged_chars, merged_chars;

  while (!is_path_start(row, col)) {
    flat_id = flat_index(col, row, ncol);
    switch(path[flat_id]) {
      case int('d'):
        col--; row--;
        unmerged_chars = get_current_chars(row, col, sequences, qualities);
        merged_chars = merge_forward_and_reverse_chars(
          unmerged_chars, merged_qualities_match, merged_qualities_mismatch
        );
        aligned.sequence.push_back(merged_chars["sequence"]);
        aligned.quality.push_back(merged_chars["quality"]);
        aligned.overlap_length++;
        if (is_same_sequence_char(unmerged_chars)) aligned.overlap_matches++;
        break;
      case int('l'):
        col--;
        aligned.sequence.push_back((sequences.forward)[col]);
        aligned.quality.push_back((qualities.forward)[col]);
        if (row > 0) aligned.overlap_length++;
        break;
      case int('u'):
        row--;
        aligned.sequence.push_back((sequences.reverse)[row]);
        aligned.quality.push_back((qualities.reverse)[row]);
        if (col < ncol - 1) aligned.overlap_length++;
        break;
      default:
        Rcpp::stop("Invalid backtracking value.");
    }
  }
  std::reverse(aligned.sequence.begin(), aligned.sequence.end());
  std::reverse(aligned.quality.begin(), aligned.quality.end());
  return aligned;
}

// [[Rcpp::export]]
MergedAlignment align_seqs_and_quals(
    Rcpp::List &sequences,
    Rcpp::List &qualities,
    Rcpp::IntegerVector &alignment_scores,
    Vector2d<unsigned int> &merged_qualities_match,
    Vector2d<unsigned int> &merged_qualities_mismatch
) {
  MergedAlignment alignment;
  Vector2d<int> scoring_matrix;
  Hashmap<std::vector<int>> score_and_path;

  PairedString sequences_pair(sequences["forward"], sequences["reverse"]);
  PairedString qualities_pair(qualities["forward"], qualities["reverse"]);
  scoring_matrix = create_scoring_matrix(
    alignment_scores["match"], alignment_scores["mismatch"]
  );
  score_and_path = find_best_scoring_path(
    sequences_pair, scoring_matrix, alignment_scores["gap_penalty"]
  );
  alignment = merge_by_path_backtrack(
    score_and_path["path"],
    sequences_pair,
    qualities_pair,
    merged_qualities_match,
    merged_qualities_mismatch
  );
  return(alignment);
}
