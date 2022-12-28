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
  return std::unordered_map<std::string, unsigned int> {
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
  // for (unsigned int i = 0; i < directions.size(); i++) {
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
  std::vector<std::string> &sequences,
  std::vector<std::vector<int>> &scoring_matrix,
  unsigned int row,
  unsigned int column,
  unsigned int ncol,
  int gap_p
) {
  char char1 = sequences[0][column];
  char char2 = sequences[1][row];
  return std::unordered_map<std::string, int> {
    {"left", gap_p},
    {"upper", (column == ncol - 1) ? 0 : gap_p},
    {"diag",  scoring_matrix[to_int(char1)][to_int(char2)]}
  };
}

void calc_score_path_other(
    std::vector<int> &score,
    std::vector<int> &path,
    std::vector<std::string> &sequences,
    std::vector<std::vector<int>> &scoring_matrix,
    int gap_p
) {
  std::unordered_map<std::string, unsigned int> index;
  std::unordered_map<std::string, int> score_prev, score_delta, score_current;
  std::string highest_scoring_dir;
  unsigned int ncol = sequences[0].size() + 1;
  unsigned int nrow = sequences[1].size() + 1;

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

// [[Rcpp::export]]
std::vector<std::unordered_map<std::string, char>> create_alignment(
    std::vector<int> &path,
    std::vector<std::string> &sequences,
    std::vector<std::string> &qualities
) {
  unsigned int ncol = sequences[0].size() + 1, nrow = sequences[1].size() + 1;
  unsigned int col = ncol - 1, row = nrow - 1, flat_id;
  std::vector<std::unordered_map<std::string, char>> aligned;
  std::unordered_map<std::string, char> chars;
  
  while (col > 0 || row > 0) {
    flat_id = flat_index(col, row, ncol);
    chars = {{"seq1", '-'}, {"seq2", '-'}, {"qua1", ' '}, {"qua2", ' '}};
    if (path[flat_id] == int('d') || path[flat_id] == int('l')) {
      col--;
      chars["seq1"] = sequences[0][col];
      chars["qua1"] = qualities[0][col];
    }
    if (path[flat_id] == int('d') || path[flat_id] == int('u')) {
      row--;
      chars["seq2"] = sequences[1][row];
      chars["qua2"] = qualities[1][row];
    }
    aligned.push_back(chars);
  }
  
  return aligned;
}

// [[Rcpp::export]]
std::unordered_map<std::string, std::vector<int>> find_best_scoring_overlap(
    std::vector<std::string> &sequences,
    std::vector<std::vector<int>> scoring_matrix,
    unsigned int gap_p
) {
  unsigned int ncol = sequences[0].size() + 1;
  unsigned int nrow = sequences[1].size() + 1;
  std::vector<int> score(nrow * ncol);
  std::vector<int> path(nrow * ncol);
  calc_score_path_first_row(score, path, ncol);
  calc_score_path_first_column(score, path, nrow, ncol, gap_p);
  calc_score_path_other(score, path, sequences, scoring_matrix, gap_p);
  return std::unordered_map<std::string, std::vector<int>> {
    {"score", score}, {"path", path}
  };
}

// # [[Rcpp::export]]
std::unordered_map<std::string, std::string> merge_alignment(
  std::vector<std::unordered_map<std::string, char>> &alignment,
  std::unordered_map<std::string, std::vector<std::vector<unsigned int>>>
    &qual_merge_map
) {
  std::unordered_map<std::string, std::string> seqs_and_quals;
  // iterate alignment in reverse
  //  check q score, then select the letter with larger score
  //  calculate posterior q
  return seqs_and_quals;
}

// #[[Rcpp::export]]
std::vector<std::unordered_map<std::string, char>> align_seqs_and_quals(
    std::vector<std::string> &sequences,
    std::vector<std::string> &qualities,
    std::unordered_map<std::string, int> &alignment_scores,
    std::vector<std::vector<unsigned int>> &qual_merge_map
    // int match,
    // int mismatch,
    // int gap_p
) {
  std::vector<std::unordered_map<std::string, char>> alignment;
  std::vector<std::vector<int>> scoring_matrix;
  std::unordered_map<std::string, std::vector<int>> score_and_path;
  
  scoring_matrix = create_scoring_matrix(
    alignment_scores["match"], alignment_scores["mismatch"]
  );
  score_and_path = find_best_scoring_overlap(
    sequences, scoring_matrix, alignment_scores["gap_penalty"]);
  alignment = create_alignment(score_and_path["path"], sequences, qualities);
  // merged_alignment = merge_alignment(alignment)
  return(alignment);
}
