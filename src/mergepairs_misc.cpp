#include <vector>

unsigned int flatten_index(
    unsigned int row,
    unsigned int column,
    unsigned int num_columns
) {
  return row * num_columns + column;
}

// [[Rcpp::export]]
int to_int(char &nt) {
  int int_nt;
  switch(nt) {
  case 'A':
    int_nt = 0;
    break;
  case 'C':
    int_nt = 1;
    break;
  case 'G':
    int_nt = 2;
    break;
  case 'T':
    int_nt = 3;
    break;
  default:
    int_nt = 4;
  }
  return(int_nt);
}
