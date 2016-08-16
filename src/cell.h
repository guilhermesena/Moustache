#include <vector>
#include <iostream>
#include <string>

using std::vector;
using std::string;

// Class that will store all useful data about the cell:
// nonzero coordinates, cluster index, name and index
class Cell {
private:
	struct GeneCoord {
		size_t gene_index;
		double val;

		GeneCoord(size_t _gene_index, double _val) :
				gene_index(_gene_index), val(_val) {
		}

		bool operator <(const GeneCoord &rhs) const{
			return gene_index < rhs.gene_index;
		}
	};

public:
	//Graph edges
	vector<Cell*> adj;
	vector<double> adj_dist;

	//Gene values
	vector<size_t> ind;
	vector<double> val;

	//Cell identifiers
	int index;
	string name;

	Cell(string, size_t);
	void InsertGene(int, double);
	void InsertEdge(Cell*, double);
	void AdjustGeneReads();
	void print();
};
