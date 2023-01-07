#include <string>

class MergedAlignment {
private:
  std::string sequence = "";
  std::string quality = "";
  unsigned int overlap_length = 0;
  unsigned int overlap_matches = 0;
public:
  MergedAlignment() {}
  std::string get_sequence() {
    return sequence;
  }
  std::string get_quality() {
    return quality;
  }
  MergedAlignment &append()
};

