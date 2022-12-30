#include <Rcpp.h>
#include <fstream>
using namespace Rcpp;

void calc_score_path_other(
    std::vector<int> &score,
    std::vector<int> &path,
    std::unordered_map<std::string, std::string> &sequences,
    std::vector<std::vector<int>> &scoring_matrix,
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

std::unordered_map<std::string, std::vector<int>> find_best_scoring_path(
    std::unordered_map<std::string, std::string> &sequences,
    std::vector<std::vector<int>> scoring_matrix,
    int gap_p
);

std::unordered_map<std::string, std::string> merge_by_path_backtrack(
    std::vector<int> &path,
    std::unordered_map<std::string, std::string> &sequences,
    std::unordered_map<std::string, std::string> &qualities,
    std::vector<std::vector<unsigned int>> &merged_qualities_match,
    std::vector<std::vector<unsigned int>> &merged_qualities_mismatch
);

std::unordered_map<std::string, char> merge_forward_and_reverse_chars(
    std::unordered_map<std::string, char> chars,
    std::vector<std::vector<unsigned int>> merged_qualities_match,
    std::vector<std::vector<unsigned int>> merged_qualities_mismatch
);