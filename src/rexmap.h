#include <Rcpp.h>
#include <fstream>
using namespace Rcpp;

// Constants
#define PHRED_OFFSET 33
#define GAP_P -7

void nt2int(char *oseq, const char *iseq);
void int2nt(char *oseq, const char *iseq);
