#include <Rcpp.h>
#include <fstream>
using namespace Rcpp;

void calc_score_path_other(
    std::vector<int> &score,
    std::vector<int> &path,
    std::vector<std::string> &sequences,
    std::vector<std::vector<int>> &align_scores,
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