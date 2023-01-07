#include <vector>
#include <unordered_map>
#include "Index.h"
#include "PairedString.h"
#include "NodeScore.h"

int to_int(char &nt);

/*
 * Class for solving sequence alignment using Needleman-Wunsch algorithm
 */
class AlignmentGrid {

private:
  
  std::string alignment_type;
  std::vector<int> score;
  std::vector<int> path;
  unsigned int num_rows;
  unsigned int num_columns;
  
  NodeScore find_previous_scores(
    unsigned int &row,
    unsigned int &column
  );
  char find_highest_scoring_direction(
    NodeScore node_scores
  );
  NodeScore calculate_score_deltas(
    PairedString &sequences,
    std::vector<std::vector<int>> &scoring_matrix,
    unsigned int &row,
    unsigned int &column,
    int &gap_penalty
  );
  NodeScore calculate_current_scores(
    unsigned int &row,
    unsigned int &column,
    PairedString &sequences,
    std::vector<std::vector<int>> &scoring_matrix,
    int &gap_penalty
  );
  
public:
  
  AlignmentGrid(
    unsigned int num_rows,
    unsigned int num_columns,
    std::string type
  );
  AlignmentGrid &fill_first_row();
  AlignmentGrid &fill_first_column(
    int &gap_penalty
  );
  AlignmentGrid &fill_other_rows_cols(
    PairedString &sequences,
    std::vector<std::vector<int>> &scoring_matrix,
    int &gap_penalty
  );
  AlignmentGrid &find_best_scoring_path(
      PairedString &sequences,
      std::vector<std::vector<int>> &scoring_matrix,
      int &gap_penalty
  );
  AlignmentGrid &set_manual_score(
    std::vector<int> score
  );
  AlignmentGrid &set_manual_path(
    std::vector<int> path
  );
  std::vector<int> get_score();
  std::vector<int> get_path();
  
};
