#include "Index.h"

Index::Index(
  unsigned int current_id,
  unsigned int left_id,
  unsigned int upper_id,
  unsigned int diag_id
) : current(current_id), left(left_id), upper(upper_id), diag(diag_id) {}

Index::Index(
  unsigned int row,
  unsigned int column,
  unsigned int num_columns
) {
  current = flatten_index(row, column, num_columns);
  left = flatten_index(row, column - 1, num_columns);
  upper = flatten_index(row - 1, column, num_columns);
  diag = flatten_index(row - 1, column - 1, num_columns);  
}