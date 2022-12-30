#include "mergepairs.h"
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

void calc_score_path_first_column(
    std::vector<int> &score,
    std::vector<int> &path,
    int nrow,
    int ncol,
    int gap_p
) {
  unsigned int flat_id;
  for (int row = 1; row < nrow; row++) {
    flat_id = flat_index(0, row, ncol);
    score[flat_id] = gap_p * row;
    path[flat_id] = int('u');
  }
  return;
}

void calc_score_path_first_row(
    std::vector<int> &score,
    std::vector<int> &path,
    int ncol
) {
  for (int flat_id = 0; flat_id < ncol; flat_id++) {
    score[flat_id] = 0;
    path[flat_id] = int('l');
  }
  return;
}

// [[Rcpp::export]]
std::unordered_map<std::string, unsigned int> get_indexes(
    unsigned int &column,
    unsigned int &row,
    unsigned int &ncol
) {
  return {
    {"current", flat_index(column, row, ncol)},
    {"left", flat_index(column - 1, row, ncol)},
    {"upper", flat_index(column, row - 1, ncol)},
    {"diag", flat_index(column - 1, row - 1, ncol)}
  };
}

std::unordered_map<std::string, int> get_left_upper_diag_scores(
    std::vector<int> &score,
    std::unordered_map<std::string, unsigned int> index
) {
  std::unordered_map<std::string, int> score_prev;
  std::vector<std::string> directions = {"left", "upper", "diag"};
  for (const std::string &direction : directions) {
    score_prev[direction] = score[index[direction]];
  }
  return score_prev;
}

// Find direction of the highest scoring in backtracking the alignment graph
//  - brake ties by preferring upper and left over diag (random choice)
std::string get_highest_scoring_direction(
  std::unordered_map<std::string, int> score_from
) {
  std::string direction;
  if (score_from["upper"] >= score_from["diag"] &&
      score_from["upper"] >= score_from["left"]) {
    direction = "upper";
  } else if (score_from["left"] >= score_from["diag"]) {
    direction = "left";
  } else {
    direction = "diag";
  }
  return direction;
}

std::unordered_map<std::string, int> get_score_deltas(
  std::unordered_map<std::string, std::string> &sequences,
  std::vector<std::vector<int>> &scoring_matrix,
  unsigned int row,
  unsigned int column,
  unsigned int ncol,
  int gap_p
) {
  char char1 = sequences["forward"][column];
  char char2 = sequences["reverse"][row];
  return {
    {"left", gap_p},
    {"upper", (column == ncol - 1) ? 0 : gap_p},
    {"diag",  scoring_matrix[to_int(char1)][to_int(char2)]}
  };
}

void calc_score_path_other(
    std::vector<int> &score,
    std::vector<int> &path,
    std::unordered_map<std::string, std::string> &sequences,
    std::vector<std::vector<int>> &scoring_matrix,
    int gap_p
) {
  std::unordered_map<std::string, unsigned int> index;
  std::unordered_map<std::string, int> score_prev, score_delta, score_current;
  std::string highest_scoring_dir;
  unsigned int ncol = sequences["forward"].size() + 1;
  unsigned int nrow = sequences["reverse"].size() + 1;

  for (unsigned int row = 1; row < nrow; row++) {
    for (unsigned int column = 1; column < ncol; column++) {
      index = get_indexes(column, row, ncol);
      score_prev = get_left_upper_diag_scores(score, index);
      score_delta = get_score_deltas(
        sequences, scoring_matrix, row, column, ncol, gap_p
      );
      score_current = sum_values(score_prev, score_delta);
      highest_scoring_dir = get_highest_scoring_direction(score_current);
      
      score[index["current"]] = score_current[highest_scoring_dir];
      path[index["current"]] = int(highest_scoring_dir.front());
    }
  }
  return;
}

std::unordered_map<std::string, char> get_current_chars(
    unsigned int row,
    unsigned int col,
    std::unordered_map<std::string, std::string> &sequences,
    std::unordered_map<std::string, std::string> &qualities
) {
  return {
    {"sequence_forward", sequences["forward"][col]},
    {"quality_forward", qualities["forward"][col]},
    {"sequence_reverse", sequences["reverse"][row]},
    {"quality_reverse", qualities["reverse"][row]}
  };
}

// [[Rcpp::export]]
char get_merged_qualities(
    char &q1char,
    char &q2char,
    std::vector<std::vector<unsigned int>> &merged_qualities
) {
  unsigned int q1 = int(q1char) - PHRED_OFFSET;
  unsigned int q2 = int(q2char) - PHRED_OFFSET;
  unsigned int q_merged;
  if (q2 > q1) std::swap(q1, q2);
  q_merged = merged_qualities[q1 - 1][q2 - 1] + PHRED_OFFSET;
  return char(q_merged);
}

std::unordered_map<std::string, char> merge_forward_and_reverse_chars(
    std::unordered_map<std::string, char> chars,
    std::vector<std::vector<unsigned int>> merged_qualities_match,
    std::vector<std::vector<unsigned int>> merged_qualities_mismatch
) {
  std::unordered_map<std::string, char> merged;
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

std::unordered_map<std::string, std::string> merge_by_path_backtrack(
    std::vector<int> &path,
    std::unordered_map<std::string, std::string> &sequences,
    std::unordered_map<std::string, std::string> &qualities,
    std::vector<std::vector<unsigned int>> &merged_qualities_match,
    std::vector<std::vector<unsigned int>> &merged_qualities_mismatch
) {
  unsigned int ncol = sequences["forward"].size() + 1;
  unsigned int nrow = sequences["reverse"].size() + 1;
  unsigned int col = ncol - 1;
  unsigned int row = nrow - 1;
  unsigned int flat_id;
  std::unordered_map<std::string, std::string> aligned;
  std::unordered_map<std::string, char> unmerged_chars, merged_chars;
  std::string merged_sequence, merged_quality;
  
  while (col > 0 || row > 0) {
    flat_id = flat_index(col, row, ncol);
    switch(path[flat_id]) {
      case int('d'):
        col--; row--;
        unmerged_chars = get_current_chars(row, col, sequences, qualities);
        merged_chars = merge_forward_and_reverse_chars(
          unmerged_chars, merged_qualities_match, merged_qualities_mismatch
        );
        aligned["sequence"].push_back(merged_chars["sequence"]);
        aligned["quality"].push_back(merged_chars["quality"]);
        break;
      case int('l'):
        col--;
        aligned["sequence"].push_back(sequences["forward"][col]);
        aligned["quality"].push_back(qualities["forward"][col]);
        break;
      case int('u'):
        row--;
        aligned["sequence"].push_back(sequences["reverse"][row]);
        aligned["quality"].push_back(qualities["reverse"][row]);
        break;
      default:
        Rcpp::stop("Invalid backtracking value.");
    }
  }
  std::reverse(aligned["sequence"].begin(), aligned["sequence"].end());
  std::reverse(aligned["quality"].begin(), aligned["quality"].end());
  return aligned;
}

std::unordered_map<std::string, std::vector<int>> find_best_scoring_path(
    std::unordered_map<std::string, std::string> &sequences,
    std::vector<std::vector<int>> scoring_matrix,
    int gap_p
) {
  unsigned int ncol = sequences["forward"].size() + 1;
  unsigned int nrow = sequences["reverse"].size() + 1;
  std::vector<int> score(nrow * ncol);
  std::vector<int> path(nrow * ncol);
  calc_score_path_first_row(score, path, ncol);
  calc_score_path_first_column(score, path, nrow, ncol, gap_p);
  calc_score_path_other(score, path, sequences, scoring_matrix, gap_p);
  return {{"score", score}, {"path", path}};
}

std::unordered_map<std::string, std::string> align_seqs_and_quals(
    std::unordered_map<std::string, std::string> &sequences,
    std::unordered_map<std::string, std::string> &qualities,
    std::unordered_map<std::string, int> &alignment_scores,
    std::vector<std::vector<unsigned int>> &merged_qualities_match,
    std::vector<std::vector<unsigned int>> &merged_qualities_mismatch
) {
  std::unordered_map<std::string, std::string> alignment;
  std::vector<std::vector<int>> scoring_matrix;
  std::unordered_map<std::string, std::vector<int>> score_and_path;
  
  scoring_matrix = create_scoring_matrix(
    alignment_scores["match"], alignment_scores["mismatch"]
  );
  score_and_path = find_best_scoring_path(
    sequences, scoring_matrix, alignment_scores["gap_penalty"]
  );
  alignment = merge_by_path_backtrack(
    score_and_path["path"],
    sequences,
    qualities,
    merged_qualities_match,
    merged_qualities_mismatch
  );
  return(alignment);
}
