#include "rexmap.h"
#include <cstdio>
// [[Rcpp::plugins(cpp11)]]

unsigned int flat_index(
    unsigned int column_index,
    unsigned int row_index,
    unsigned int ncol) {
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
  for (int row_index = 1; row_index < nrow; row_index++) {
    flat_id = flat_index(0, row_index, ncol);
    score[flat_id] = gap_p * row_index;
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

std::unordered_map<std::string, int> sum(
    std::unordered_map<std::string, int> map1,
    std::unordered_map<std::string, int> map2
) {
  std::unordered_map<std::string, int> sum;
  
  for (const auto &key_value : map1) {
    sum[key_value.first] = key_value.second;
  }
  
  for (const auto &key_value : map2) {
    sum[key_value.first] += key_value.second;
  }
  
  return sum;
}

std::unordered_map<std::string, int> get_left_upper_diag_scores(
    std::vector<int> &score,
    std::unordered_map<std::string, unsigned int> &index
) {
  std::unordered_map<std::string, int> score_prev;
  std::vector<std::string> directions = {"left", "upper", "diag"};
  for (const std::string &direction : directions) {
    score_prev[direction] = index[direction];
  }
  return score_prev;
}

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
  std::string &seq1,
  std::string &seq2,
  std::vector<std::vector<int>> &align_scores,
  unsigned int row,
  unsigned int column,
  unsigned int ncol,
  int gap_p
) {
  return std::unordered_map<std::string, int> {
    {"left", gap_p},
    {"upper", (column == ncol - 1) ? 0 : gap_p},
    {"diag",  align_scores[to_int(seq1[column])][to_int(seq2[row])]}
  };
}

void calc_score_path_other(
    std::vector<int> &score,
    std::vector<int> &path,
    std::string &seq1,
    std::string &seq2,
    std::vector<std::vector<int>> &align_scores,
    unsigned int nrow,
    unsigned int ncol,
    int gap_p
) {
  std::unordered_map<std::string, unsigned int> index;
  std::unordered_map<std::string, int> score_prev, score_delta, score_current;
  std::string highest_scoring_dir;

  for (unsigned int row = 1; row < nrow; row++) {
    for (unsigned int column = 1; column < ncol; column++) {
      index = get_indexes(column, row, ncol);
      score_prev = get_left_upper_diag_scores(score, index);
      score_delta = get_score_deltas(
        seq1, seq2, align_scores, row, column, ncol, gap_p
      );
      score_current = sum(score_prev, score_delta);
      highest_scoring_dir = get_highest_scoring_direction(score_current);
      
      score[index["current"]] = score_current[highest_scoring_dir];
      path[index["current"]] = int(highest_scoring_dir.front());
    }
  }
  return;
}

// #[[Rcpp::export]]
std::vector<std::unordered_map<std::string, char>> create_alignment_from_path(
    std::vector<char> &path,
    std::string &seq1,
    std::string &seq2,
    std::string &qua1,
    std::string &qua2,
    unsigned int nrow,
    unsigned int ncol
) {
  std::vector<std::unordered_map<std::string, char>> aligned;
  // Trace back over p to form the alignment.
  size_t aligned_index = 0;
  unsigned int column_index = ncol - 1;
  unsigned int row_index = nrow - 1;
  
  while (column_index > 0 || row_index > 0) {
    switch (path[column_index * ncol + row_index]) {
    case 'd':
      column_index--;
      row_index--;
      aligned[aligned_index] = {
        {"seq1", seq1[column_index]},
        {"seq2", seq2[row_index]},
        {"qua1", qua1[column_index]},
        {"qua2", qua2[row_index]}
      };
      break;
    case 'u':
      row_index--;
      aligned[aligned_index] = {
        {"seq1", '-'},
        {"seq2", seq2[row_index]},
        {"qua1", ' '},
        {"qua2", qua2[row_index]}
      };
      break;
    case 'l':
      column_index--;
      aligned[aligned_index] = {
        {"seq1", seq1[column_index]},
        {"seq2", '-'},
        {"qua1", qua1[column_index]},
        {"qua2", ' '}
      };
      break;
    default:
      Rcpp::stop("N-W Align out of range.");
    }
    aligned_index++;
  }
  return aligned;
}

// [[Rcpp::export]]
std::unordered_map<std::string, std::vector<int>> find_overlap_alignment(
    std::string &seq1,
    std::string &seq2,
    unsigned int nrow,
    unsigned int ncol,
    std::vector<std::vector<int>> scoring_matrix,
    unsigned int gap_p
) {
  std::vector<int> score(nrow * ncol);
  std::vector<int> path(nrow * ncol);
  calc_score_path_first_row(score, path, ncol);
  calc_score_path_first_column(score, path, nrow, ncol, gap_p);
  calc_score_path_other(
    score, path, seq1, seq2, scoring_matrix, nrow, ncol, gap_p
  );
  
  return std::unordered_map<std::string, std::vector<int>> {
    {"score", score}, {"path", path}
  };
}

// [[Rcpp::export]]
std::vector<std::unordered_map<std::string, char>> align_seqs_and_quals(
    std::string &seq1,
    std::string &seq2,
    std::string &qua1,
    std::string &qua2,
    int match,
    int mismatch,
    int gap_p
) {
  std::vector<std::unordered_map<std::string, char>> alignment;

  unsigned int ncol = seq1.size() + 1;
  unsigned int nrow = seq2.size() + 1;
  std::vector<std::vector<int>> scoring_matrix = create_scores(match, mismatch);
  Rcpp::List score_and_path = find_overlap_alignment(
    seq1, seq2, nrow, ncol, scoring_matrix, gap_p
  );
  alignment = create_alignment_from_path(
    score_and_path["path"], seq1, seq2, qua1, qua2, nrow, ncol
  );
  // merged_alignment = merge_alignment(alignment)
  return(alignment);
}

// [[Rcpp::export]]
void test_fun(std::vector<int> x) {
  for (unsigned int i = x.size(); i-- > 0; ) {
    std::cout << "x[i] = " << x[i] << "\n";
  }
  return;
}
