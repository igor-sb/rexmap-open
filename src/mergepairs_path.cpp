#include <Rcpp.h>
#include "rexmap.h"
// 
// void calc_score_path_first_column(
//     std::vector<int> &score,
//     std::vector<int> &path,
//     int nrow,
//     int ncol,
//     int gap_p
// ) {
//   unsigned int flat_id;
//   for (int row = 1; row < nrow; row++) {
//     flat_id = flatten_index(0, row, ncol);
//     score[flat_id] = gap_p * row;
//     path[flat_id] = int('u');
//   }
//   return;
// }
// 
// void calc_score_path_first_row(
//     std::vector<int> &score,
//     std::vector<int> &path,
//     int ncol
// ) {
//   for (int flat_id = 0; flat_id < ncol; flat_id++) {
//     score[flat_id] = 0;
//     path[flat_id] = int('l');
//   }
//   return;
// }
// 
// NodeScore find_previous_scores(
//     std::vector<int> &score,
//     Index index
// ) {
//   NodeScore previous_scores(
//     score[index.left], score[index.upper], score[index.diag]
//   );
//   return previous_scores;
// }
// 
// // Find direction of the highest scoring in backtracking the alignment graph
// //  - brake ties by preferring upper and left over diag (random choice)
// char find_highest_scoring_direction(
//     NodeScore previous_scores
// ) {
//   char direction;
//   if (previous_scores.upper >= previous_scores.diag &&
//       previous_scores.upper >= previous_scores.left) {
//     return 'u';
//   } else if (previous_scores.left >= previous_scores.diag) {
//     return 'l';
//   } else {
//     return 'd';
//   }
// }
// 
// NodeScore calculate_score_deltas(
//     PairedString &sequences,
//     Vector2d<int> &scoring_matrix,
//     unsigned int row,
//     unsigned int column,
//     unsigned int ncol,
//     int gap_p
// ) {
//   char char1 = (sequences.forward)[column - 1];
//   char char2 = (sequences.reverse)[row - 1];
//   NodeScore score_deltas;
//   score_deltas.left = gap_p;
//   score_deltas.upper = (column == ncol - 1) ? 0 : gap_p;
//   score_deltas.diag = scoring_matrix[to_int(char1)][to_int(char2)];
//   return score_deltas;
// }
// 
// void calc_score_path_other(
//     std::vector<int> &score,
//     std::vector<int> &path,
//     PairedString &sequences,
//     Vector2d<int> &scoring_matrix,
//     int &gap_p
// ) {
//   Index index;
//   NodeScore scores_previous, scores_delta, scores_current;
//   char highest_scoring_dir;
//   unsigned int ncol = sequences.forward.size() + 1;
//   unsigned int nrow = sequences.reverse.size() + 1;
//   
//   for (unsigned int row = 1; row < nrow; row++) {
//     for (unsigned int column = 1; column < ncol; column++) {
//       index = create_flat_indexes(column, row, ncol);
//       scores_previous = find_previous_scores(score, index);
//       scores_delta = calculate_score_deltas(
//         sequences, scoring_matrix, row, column, ncol, gap_p
//       );
//       scores_current = add_scores(scores_previous, scores_delta);
//       highest_scoring_dir = find_highest_scoring_direction(scores_current);
// 
//       score[index.current] = scores_current.get_direction(highest_scoring_dir);
//       path[index.current] = int(highest_scoring_dir);
//     }
//   }
//   return;
// }
// 
// Hashmap<std::vector<int>> find_best_scoring_path(
//     PairedString &sequences,
//     Vector2d<int> scoring_matrix,
//     int gap_p
// ) {
//   unsigned int ncol = sequences.forward.size() + 1;
//   unsigned int nrow = sequences.reverse.size() + 1;
//   std::vector<int> score(ncol * nrow), path(ncol * nrow);
//   calc_score_path_first_row(score, path, ncol);
//   calc_score_path_first_column(score, path, nrow, ncol, gap_p);
//   calc_score_path_other(score, path, sequences, scoring_matrix, gap_p);
//   return {{"score", score}, {"path", path}};
// }
// 
// AlignmentGrid find_best_scoring_path2(
//     PairedString &sequences,
//     Vector2d<int> scoring_matrix,
//     int gap_penalty
// ) {
//   unsigned int ncol = sequences.forward.size() + 1;
//   unsigned int nrow = sequences.reverse.size() + 1;
//   AlignmentGrid score_and_path(nrow, ncol);
//   score_and_path
//     .fill_first_row()
//     .fill_first_column(gap_penalty)
//     .fill_other_rows_cols(sequences, scoring_matrix, gap_penalty);
//   return score_and_path;
// }