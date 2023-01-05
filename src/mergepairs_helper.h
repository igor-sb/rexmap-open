#include <Rcpp.h>
#include "rexmap.h"

void calc_score_path_other(
    std::vector<int> &score,
    std::vector<int> &path,
    PairedString &sequences,
    Vector2d<int> &scoring_matrix,
    int gap_p
);

void calc_score_path_first_row(
    std::vector<int> &score,
    std::vector<int> &path,
    int ncol
);

void calc_score_path_first_column(
    std::vector<int> &score,
    std::vector<int> &path,
    int nrow,
    int ncol,
    int gap_p
);

Hashmap<std::vector<int>> find_best_scoring_path(
    PairedString &sequences,
    Vector2d<int> scoring_matrix,
    int gap_p
);

MergedAlignment merge_by_path_backtrack(
    std::vector<int> &path,
    PairedString &sequences,
    PairedString &qualities,
    Vector2d<unsigned int> &merged_qualities_match,
    Vector2d<unsigned int> &merged_qualities_mismatch
);

Hashmap<char> merge_forward_and_reverse_chars(
    Hashmap<char> chars,
    Vector2d<unsigned int> merged_qualities_match,
    Vector2d<unsigned int> merged_qualities_mismatch
);