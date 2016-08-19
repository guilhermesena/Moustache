#include "cellset.h"
#include <sstream>
#include <string>
#include <iostream>
#include <utility>
#include <algorithm>
#define eps 1e-6

using std::cerr;
using std::cout;
using std::string;
using std::istringstream;
using std::min;
using std::pair;

CellSet::CellSet(size_t _NumNeighbors, bool _ExcludeNoise, istream& _st,
		ostream &_out) :
		NumNeighbors(_NumNeighbors), ExcludeNoise(_ExcludeNoise), st(_st), out(
				_out) {

	dist_memory.clear();
	numgenes = 0;
}

void CellSet::FreeCells() {
	for (int i = 0; i < cells.size(); i++) {
		delete cells[i];
	}
}

CellSet::~CellSet() {
	FreeCells();
}

unordered_map<pair<size_t, size_t>, double> CellSet::dist_memory;

//Reads the drop-seq dataset from an istream chosen by the user (ie, stdin)
void CellSet::ReadFromStream() {
	string cellname, line, genename;
	float count;
	size_t ncells = 0;

	getline(st, line);
	istringstream rownames(line);

	while (rownames >> cellname) {
		cells.push_back(new Cell(cellname, (int) ncells++));
	}

	while (st >> genename) {
		for (int i = 0; i < ncells; i++) {
			st >> count;
			if (count > eps)
				cells[i]->InsertGene(numgenes, count);
		}
		genenames[numgenes++] = genename;
	}

	//If number of neighbors is bigger than number of cells, update it
	NumNeighbors = min(NumNeighbors, cells.size() - 1);

#ifdef DEBUG
	cerr << "total cells: " << cells.size() << "\n";
	cerr << "total genes: " << numgenes << "\n";
#endif

}

void CellSet::CreateVpTree() {
	vpt.create(&cells);

#ifdef DEBUG
	PrintCells();
#endif

}

void CellSet::BuildNNGraph() {
	BuildNNGraph(ExcludeNoise);
}

//Builds edges to k nearest neighbors
void CellSet::BuildNNGraph(bool _ExcludeNoise) {

	vector<NNPair> dists;
	size_t ncells = cells.size();
	vector<pair<vector<Cell*>, vector<double> > > list_edges;
	for (size_t i = 0; i < ncells; i++) {
		vector<Cell*> res;
		vector<double> dst;
		vpt.search(cells[i], NumNeighbors + 1, &res, &dst);

		list_edges.push_back(make_pair(res, dst));
		if (_ExcludeNoise) {
			dists.push_back(NNPair(i, dst[min(NumNeighbors, (size_t) 4)]));
		}
	}

	//If we don't need to exclude noise cells, we're done
	if (!_ExcludeNoise)
		return;

	cerr << "Removing outlier cells...\n";
	RemoveOutlierCells(dists);

	//stores remaining cells
	unordered_map<Cell*, int> reg;
	for (auto c : cells) {
		reg[c] = 1;
	}

	//Adds edges that do not contain outlier cells
	for (size_t i = 0; i < cells.size(); i++) {
		for (size_t j = 0; j <= NumNeighbors; j++) {
			Cell *cc = list_edges[i].first[j];

			//Add edge if cell has not been excluded
			if (reg[cc]) {
				cells[i]->InsertEdge(cc, list_edges[i].second[j]);
			}
		}
	}
}

//Removes cells whose k-th nearest neighbor is too far
void CellSet::RemoveOutlierCells(vector<NNPair> &dists) {
	sort(dists.begin(), dists.end());

	double mn = dists[0].distance, mx = dists[dists.size() - 1].distance;

	//Removes obvious outliers from count
	for (size_t j = dists.size() - 1;
			j > 0 && dists[j].distance >= 10.0 * dists[j - 1].distance; j--)
		mx = dists[j - 1].distance;

	double threshold = -1.0;
	size_t oldcellsize = cells.size(), ind_threshold;

	//Assumption: At least 80% of cells are not noise. Look for spike in other 20%
	for (size_t j = 4* oldcellsize/5; j < oldcellsize; j++) {
		if (dists[j].distance - dists[j - 1].distance
				>= 2.0 * (mx - mn) / ((double) oldcellsize)) {
			threshold = dists[j - 1].distance;
			ind_threshold = j - 1;

			//Discard noise cells
			vector<size_t> discard;
			for (size_t k = j + 1; k < oldcellsize; k++) {
				discard.push_back(dists[k].CellIndex);
			}

			sort(discard.begin(), discard.end());

			//Use int bc size_t is unsigned!
			for (int k = discard.size() - 1; k >= 0; k--) {
				cells.erase(cells.begin() + discard[k]);
			}

			break;
		}
	}

	//If number of neighbors is bigger than number of cells, update it
	NumNeighbors = min(NumNeighbors, cells.size() - 1);

#ifdef DEBUG
	if (threshold == -1.0) {
		cerr << "No noise cells in the dataset\n";
	} else {
		cerr << "Distance threshold: " << threshold << "\n";
		cerr << 100.0 * (oldcellsize - ind_threshold) / ((double) oldcellsize)
		<< "\% of cells are noise\n";

		cerr << "New graph has " << cells.size() << " cells\n";
	}
#endif

}

void CellSet::MedianFilter() {
	vector<Cell*> newData;
	for (auto cell : cells) {
		Cell *nc = new Cell(cell->name, cell->index);
		*nc = cell->AdjustGeneReads();
		newData.push_back(nc);
	}

	FreeCells();
	cells = newData;
	CreateVpTree();

	//Rebuilds graph not excluding outlier cells
	BuildNNGraph(ExcludeNoise);

}

//Print some cell rows for debug purposes
void CellSet::PrintCells() {
	int i = 0;
	for (auto c : cells) {
		c->print();
		if (++i == 20)
			break;
	}
}

//Writes cell data to given ostream
void CellSet::WriteToStream() {
	//Write a row with cell names
	for (auto c : cells) {
		out << "\t" << c->name;
	}
	out << "\n";

	//Pointers to the current gene in each cell
	vector<size_t> pointers(cells.size(), 0);
	for (size_t gene = 0; gene < numgenes; gene++) {
		out << genenames[gene] << "\t";

		//Print the gene if cell has it, 0 otherwise
		for (size_t i = 0; i < cells.size(); i++) {
			if (pointers[i] < cells[i]->ind.size()
					&& cells[i]->ind[pointers[i]] == gene) {
				out << cells[i]->val[pointers[i]++] << "\t";
			} else
				out << "0\t";
		}
		out << "\n";
	}
}

void CellSet::CountConnectedComponents() {
	size_t NumConnectedComponents = 0;
	unordered_map<size_t, int> marked;

	for (size_t i = 0; i < cells.size(); i++) {
		if (!marked[cells[i]->index]) {
			NumConnectedComponents++;
			size_t CCSize = DFS(cells[i], marked);
			cerr << "CC with " << CCSize << " cells\n";
		}
	}

	cerr << "Graph has " << NumConnectedComponents << " CCs\n";
}

size_t CellSet::DFS(Cell *cell, unordered_map<size_t, int> &marked) {
	marked[cell->index] = 1;
	size_t ans = 1;
	for (auto c : cell->adj) {
		if (!marked[c->index]) {
			ans += DFS(c, marked);
		}
	}

	for (auto c : cell->from) {
		if (!marked[c->index]) {
			ans += DFS(c, marked);
		}
	}

	return ans;
}
