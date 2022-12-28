#include <Rcpp.h>
#include <fstream>
using namespace Rcpp;

// Constants
#define PHRED_OFFSET 33

void nt2int(char *oseq, const char *iseq);
void int2nt(char *oseq, const char *iseq);
int to_int(char &nt);
std::vector< std::vector<int> > create_scoring_matrix(int match, int mismatch);
std::unordered_map<std::string, int> sum_values(
    std::unordered_map<std::string, int> map1,
    std::unordered_map<std::string, int> map2
);
void print_flat_vector(const std::vector<int> &v, unsigned int ncol);