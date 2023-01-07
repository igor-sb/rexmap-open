#include <string>
#include <unordered_map>
#include <vector>
// [[Rcpp::plugins(cpp11)]]

void nt2int(char *oseq, const char *iseq);
void int2nt(char *oseq, const char *iseq);

template<typename T>
using Hashmap = std::unordered_map<std::string, T>;

#define PHRED_OFFSET 33