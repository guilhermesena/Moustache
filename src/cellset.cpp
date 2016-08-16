#include "cellset.h"
#include <sstream>
#include <string>
#include <iostream>
#define eps 1e-6

using std::cerr;
using std::cout;
using std::string;
using std::istringstream;

CellSet::CellSet() {
	dist_memory.clear();
}
CellSet::~CellSet() {
	for (int i = 0; i < cells.size(); i++) {
		delete cells[i];
	}
}

unordered_map<pair<size_t, size_t>, double> CellSet::dist_memory;

//Reads the drop-seq dataset from an istream chosen by the user (ie, stdin)
void CellSet::ReadFromStream(istream &st) {
	string cellname, line, genename;
	float count;
	size_t ncells = 0, gene = 0;

	getline(st, line);
	istringstream rownames(line);

	while (rownames >> cellname) {
		cells.push_back(new Cell(cellname, (int) ncells++));
	}

	while (st >> genename) {
		for (int i = 0; i < ncells; i++) {
			st >> count;
			if (count > eps)
				cells[i]->InsertGene(gene, count);
		}
		gene++;
	}

	cerr << "total cells: " << cells.size() << "\n";
	cerr << "total genes: " << gene << "\n";
	for (int i = 0; i < min((size_t) 20, ncells); i++)
		cells[i]->print();
}

void CellSet::CreateVpTree() {
	vpt.create(&cells);

	for (int i = 0; i < min((size_t) 20, cells.size()); i++)
		cells[i]->print();

}

void CellSet::BuildNNGraph(int NumNeighbors, bool ExcludeNoise,
		size_t *ind_thresh) {

	vector<NNPair> dists;
	size_t ncells = cells.size();

	for (size_t i = 0; i < ncells; i++) {
		vector<Cell*> res;
		vector<double> dst;
		vpt.search(cells[i], NumNeighbors + 1, &res, &dst);

		/*
		 cerr << cells[i]->name << " ";
		 for (size_t d = 0; d < dst.size(); d++)
		 cerr << dst[d] << " ";

		 cerr << "\n";
		 */

		for (size_t i = 1; i < NumNeighbors; i++) {
			cells[i]->InsertEdge(res[i], dst[i]);
		}

		if (ExcludeNoise) {
			dists.push_back(NNPair(i, dst[4]));
		}
	}

	if (ExcludeNoise) {
		sort(dists.begin(), dists.end());
		for (int i = 0; i < dists.size(); i++) {
			cout << dists[i].distance << " ";
		}

		double mn = dists[0].distance, mx = dists[dists.size() - 1].distance;

		//Removes obvious outliers from count
		for (size_t j = dists.size() - 1;
				j > 0 && dists[j].distance >= 10.0 * dists[j - 1].distance; j--)
			mx = dists[j-1].distance;

		double threshold = -1.0;
		size_t oldcellsize = cells.size();
		for (size_t j = 3 * ncells / 5; j < ncells; j++) {
			if (dists[j].distance - dists[j - 1].distance
					>= 2.0 * (mx - mn) / ((double) ncells)) {
				threshold = dists[j - 1].distance;
				if (ind_thresh != NULL)
					*ind_thresh = j - 1;

				//Discard noise cells
				vector < size_t > discard;
				for (size_t k = j + 1; k < ncells; k++) {
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
		if (threshold == -1.0) {
			cerr << "No noise cells in the dataset\n";
		} else {
			cerr << "Distance threshold: " << threshold
					<< " corresponding to index " << *ind_thresh << "\n";
			cerr
					<< 100.0 * (oldcellsize - int(*ind_thresh))
							/ ((double) oldcellsize)
					<< "\% of cells are noise\n";

			cerr << "New graph has " << cells.size() << " cells\n";
		}
	}
}

void CellSet::EstimateReads() {
	for (auto cell : cells) {
		for (auto neighbor : cell->adj) {

		}
	}
}
