class NodeScore {
public:
  int left = 0;
  int upper = 0;
  int diag = 0;
  NodeScore() {}
  NodeScore(int left_val, int upper_val, int diag_val);
  int get_at_direction(char direction);
  static NodeScore add_scores(NodeScore score1, NodeScore score2);
};