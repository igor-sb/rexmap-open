class Index {
public:
  unsigned int current;
  unsigned int left;
  unsigned int upper;
  unsigned int diag;
  Index() {}
  Index(unsigned int current_id,
        unsigned int left_id,
        unsigned int upper_id,
        unsigned int diag_id);
  Index(unsigned int row, unsigned int column, unsigned int num_columns);
};

unsigned int flatten_index(
    unsigned int row,
    unsigned int column,
    unsigned int num_columns
);
