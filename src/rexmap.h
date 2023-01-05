#include <string>
#include <unordered_map>
// [[Rcpp::plugins(cpp11)]]

class MergedAlignment {
  public:
    std::string sequence;
    std::string quality;
    unsigned int overlap_length;
    unsigned int overlap_matches;
};

class PairedString {
  public:
    std::string forward;
    std::string reverse;
    PairedString() {}
    PairedString(std::string fwd, std::string rev) {
      forward = fwd;
      reverse = rev;
    }
};

class Score {
  public:
    int left;
    int upper;
    int diag;
    Score(int left_val, int upper_val, int diag_val) :
      left(left_val), upper(upper_val), diag(diag_val) {}
};

// Aliases for common types
template<typename T>
using Hashmap = std::unordered_map<std::string, T>;

template<typename T>
using Vector2d = std::vector<std::vector<T>>;

// Constants
#define PHRED_OFFSET 33

void nt2int(char *oseq, const char *iseq);
void int2nt(char *oseq, const char *iseq);
unsigned int flat_index(
    unsigned int column_index,
    unsigned int row_index,
    unsigned int ncol
);
int to_int(char &nt);
Vector2d<int> create_scoring_matrix(int match, int mismatch);
Hashmap<int> sum_values(
    Hashmap<int> map1,
    Hashmap<int> map2
);
void print_flat_vector(const std::vector<int> &v, unsigned int ncol);
Hashmap<std::vector<int>> find_best_scoring_path(
    PairedString &sequences,
    Vector2d<int> scoring_matrix,
    int gap_p
);