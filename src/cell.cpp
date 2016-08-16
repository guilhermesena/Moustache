#include "cell.h"
#include <algorithm>
#include <set>

using std::set;
using std::cerr;
using std::min;

// Class that will store all useful data about the cell:
// nonzero coordinates, cluster index, name and index
Cell::Cell(string _name, size_t _index) :
		name(_name), index(_index) {
}

//Adds a neighbour to cell
void Cell::InsertEdge(Cell *_adj, double _adj_dist) {
	adj.push_back(_adj);
	adj_dist.push_back(_adj_dist);
}

//Adds a nonzero gene coordinate
void Cell::InsertGene(int indv, double valv) {
	ind.push_back(indv);
	val.push_back(valv);
}

//Use neighbors to adjust read values using median filter
void Cell::AdjustGeneReads() {
	size_t nn = adj.size();
	set<GeneCoord> merged_genes;

	//Add the cell's genes
	for(int i = 0; i < ind.size(); i++) {
		merged_genes.insert(GeneCoord(ind[i], val[i]));
	}

	//Add the neighbor's genes
	for (auto neighbor : adj) {
		for(int i = 0; i < neighbor->ind.size(); i++) {
			merged_genes.insert(GeneCoord(neighbor->ind[i], neighbor->val[i]));
		}
	}

	//Now we get the median gene by gene and rewrite the cell
	ind.clear();
	val.clear();
	while(!merged_genes.empty()) {
		vector<double> values;
		GeneCoord g = *(merged_genes.begin());
		merged_genes.erase(merged_genes.begin());
		values.push_back(g.val);

		//Get all the values of the same gene in all cells
		while(g.gene_index == (*(merged_genes.begin())).gene_index) {
			values.push_back((*merged_genes.begin()).val);
			merged_genes.erase(merged_genes.begin());
		}

		//Get the median if there are more than half nonzero values
		if(values.size() > nn/2) {
			sort(values.begin(), values.end());
			InsertGene(g.gene_index, values[values.size() - 1 - nn/2]);
		}
	}
}

//Prints summary of cell data (for debug purposes)
void Cell::print() {
	int sz = ind.size();
	cerr << name << "\t";
	cerr << sz << "\t";
	for (int i = 0; i < min(sz, 10); i++) {
		cerr << ind[i] << " " << val[i] << "\t";
	}
	cerr << "\n";
}
