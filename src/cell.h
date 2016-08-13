#include <vector>
#include <iostream>
#include <string>

using std::vector;
using std::string;

// Class that will store all useful data about the cell:
// nonzero coordinates, cluster index, name and index
class Cell {
public:
	vector<size_t> ind;
	vector<double> val;
	string name;
	int index;

	Cell(string _name, size_t _index);
	void insert(int indv, double valv);
	void print();
};
