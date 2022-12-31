#include "mergepairs.h"

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