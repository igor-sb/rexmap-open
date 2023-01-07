#include "NodeScore.h"
#include <stdexcept>

NodeScore::NodeScore(int left_val, int upper_val, int diag_val) :
  left(left_val), upper(upper_val), diag(diag_val) {}

int NodeScore::get_at_direction(char direction) {
  if (direction == 'l') {
    return left;
  } else if (direction == 'u') {
    return upper;
  } else if (direction == 'd') {
    return diag;
  } else {
    throw std::invalid_argument("Invalid direction");
  }
}

NodeScore NodeScore::add_scores(NodeScore scores1, NodeScore scores2) {
  return NodeScore(
    scores1.left + scores2.left,
    scores1.upper + scores2.upper,
    scores1.diag + scores2.diag
  );
}