#include <vector>
#include <iostream>
#include <string>

using std::vector;
using std::string;

// Class that will store all useful data about the cell:
// nonzero coordinates, cluster index, name and index
class Cell {
private:
public:
	//Graph edges
	vector<Cell*> adj;
	vector<double> adj_dist;

	vector<size_t> ind;
	vector<double> val;
	int index;
	string name;

	Cell(string, size_t);
	void InsertGene(int, double);
	void InsertEdge(Cell*, double);
	void print();
};
