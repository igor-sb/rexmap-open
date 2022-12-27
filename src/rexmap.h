#include <Rcpp.h>
#include <fstream>
using namespace Rcpp;

// Constants
#define PHRED_OFFSET 33

int to_int(char &nt);
std::vector< std::vector<int> > create_scores(int match, int mismatch);