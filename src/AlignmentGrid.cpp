#include "AlignmentGrid.h"

NodeScore AlignmentGrid::find_previous_scores(
    unsigned int &row,
    unsigned int &column
) {
  Index flat_ids(row, column, num_columns);
  NodeScore previous_scores(
      score[flat_ids.left], score[flat_ids.upper], score[flat_ids.diag]
  );
  return previous_scores;
}

char AlignmentGrid::find_highest_scoring_direction(
    NodeScore node_scores
) {
  if (node_scores.upper >= node_scores.diag &&
      node_scores.upper >= node_scores.left) {
    return('u');
  } else if (node_scores.left >= node_scores.diag) {
    return('l');
  } else {
    return('d');
  }
}

NodeScore AlignmentGrid::calculate_score_deltas(
    PairedString &sequences,
    std::vector<std::vector<int>> &scoring_matrix,
    unsigned int &row,
    unsigned int &column,
    int &gap_penalty
) {
  char fwd_char = (sequences.forward)[column - 1];
  char rev_char = (sequences.reverse)[row - 1];
  NodeScore score_deltas;
  score_deltas.left = gap_penalty;
  score_deltas.upper = (column == num_columns - 1) ? 0 : gap_penalty;
  score_deltas.diag = scoring_matrix[to_int(fwd_char)][to_int(rev_char)];
  return score_deltas;
}

NodeScore AlignmentGrid::calculate_current_scores(
    unsigned int &row,
    unsigned int &column,
    PairedString &sequences,
    std::vector<std::vector<int>> &scoring_matrix,
    int &gap_penalty
) {
  NodeScore scores_previous = find_previous_scores(row, column);
  NodeScore scores_delta = calculate_score_deltas(
    sequences, scoring_matrix, row, column, gap_penalty
  );
  return NodeScore::add_scores(scores_previous, scores_delta);
}

AlignmentGrid::AlignmentGrid(
    unsigned int num_rows,
    unsigned int num_columns,
    std::string type
) : score(num_rows * num_columns),
    path(num_rows * num_columns),
    num_rows(num_rows),
    num_columns(num_columns),
    alignment_type(type) {
  // check if type is "overlap" or "fitting", if not throw error
}

AlignmentGrid &AlignmentGrid::fill_first_row() {
  for (unsigned int flat_id = 0; flat_id < num_columns; flat_id++) {
    score[flat_id] = 0;
    path[flat_id] = int('l');
  }
  return *this;
}

AlignmentGrid &AlignmentGrid::fill_first_column(
    int &gap_penalty
) {
  unsigned int flat_id;
  int first_column_gap = alignment_type == "fitting" ? 0 : gap_penalty;
  for (int row = 1; row < num_rows; row++) {
    flat_id = flatten_index(row, 0, num_columns);
    score[flat_id] = first_column_gap * row;
    path[flat_id] = int('u');
  }
  return *this;
}

AlignmentGrid &AlignmentGrid::fill_other_rows_cols(
    PairedString &sequences,
    std::vector<std::vector<int>> &scoring_matrix,
    int &gap_penalty
) {
  unsigned int flat_id;
  NodeScore scores_previous, scores_delta, scores_current;
  char highest_scoring_dir;

  for (unsigned int row = 1; row < num_rows; row++) {
    for (unsigned int column = 1; column < num_columns; column++) {
      flat_id = flatten_index(row, column, num_columns);
      scores_current = calculate_current_scores(
        row, column, sequences, scoring_matrix, gap_penalty
      );
      highest_scoring_dir = find_highest_scoring_direction(scores_current);
      score[flat_id] = scores_current.get_at_direction(highest_scoring_dir);
      path[flat_id] = int(highest_scoring_dir);
    }
  }
  return *this;
}

AlignmentGrid &AlignmentGrid::fill_full_grid(
    PairedString &sequences,
    std::vector<std::vector<int>> &scoring_matrix,
    int &gap_penalty
) {
  fill_first_row();
  fill_first_column(gap_penalty);
  fill_other_rows_cols(sequences, scoring_matrix, gap_penalty);
  return *this;
}

AlignmentGrid &AlignmentGrid::set_manual_path(
    std::vector<int> path
) {
  this->path = path;
  return *this;
}

AlignmentGrid &AlignmentGrid::set_manual_score(
    std::vector<int> score
) {
  this->score = score;
  return *this;
}

std::vector<int> AlignmentGrid::get_score() {
  return score;
}

std::vector<int> AlignmentGrid::get_path() {
  return path;
}
