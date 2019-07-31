#include <vector>

class MpcPredictor
{

public:
	  MpcPredictor(double);

	  void predict(std::vector<float>&, std::vector<float>&, std::vector<float>&);

private:
	  int N = 10;
	  double T;
	  double dt;
	  //auto opti;
};

