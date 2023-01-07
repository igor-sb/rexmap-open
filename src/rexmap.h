#include <string>
#include <unordered_map>
#include <vector>
// [[Rcpp::plugins(cpp11)]]

class PairedString {
  public:
    std::string forward = "";
    std::string reverse = "";
    PairedString() {}
    PairedString(std::string fwd, std::string rev) {
      forward = fwd;
      reverse = rev;
    }
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
int to_int(char &nt);
Vector2d<int> create_scoring_matrix(int match, int mismatch);
void print_flat_vector(const std::vector<int> &v, unsigned int ncol);
Hashmap<std::vector<int>> find_best_scoring_path(
    PairedString &sequences,
    Vector2d<int> scoring_matrix,
    int gap_p
);