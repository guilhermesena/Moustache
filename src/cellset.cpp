#include "cellset.h"
#include <sstream>
#include <string>
#include <sstream>
#define eps 1e-6

using std::string;
using std::istringstream;

CellSet::CellSet() {
	dist_memory.clear();
}
CellSet::~CellSet() {
}

unordered_map<pair<size_t, size_t>, double> CellSet::dist_memory;

//Reads the drop-seq dataset from an istream chosen by the user (ie, stdin)
void CellSet::ReadFromStream(istream &st) {
	string cellname, line;
	char genename[100];
	float count;
	size_t ncells = 0, gene = 0;

	getline(st, line);
	istringstream rownames(line);

	while (rownames >> cellname) {
		cells.push_back(Cell(cellname, (int) ncells++));
	}

	while (st >> genename) {
		for (int i = 0; i < ncells; i++) {
			st >> count;
			if (count > eps)
				cells[i].insert(gene, count);
		}
		gene++;
	}

	//Debug
	cerr << "total cells: " << cells.size() << "\n";
	cerr << "total genes: " << gene << "\n";
	for (int i = 0; i < min((size_t) 20, ncells); i++)
		cells[i].print();
}

void CellSet::CreateVpTree() {
	vpt.create(&cells);

	//Debug
	for (int i = 0; i < min((size_t) 20, cells.size()); i++)
		cells[i].print();
}

double CellSet::GetDistanceThreshold(size_t *ind_thresh) {

	vector<double> dists;
	size_t ncells = cells.size();

	for (int i = 0; i < ncells; i++) {
		vector<Cell> res;
		vector<double> dst;
		vpt.search(cells[i], 5, &res, &dst);

		/*Debug
		 cerr << cells[i].name << " ";
		 for (size_t d = 0; d < dst.size(); d++) {
		 cerr << dst[d] << " ";
		 }
		 cerr << "\n";
		 */

		dists.push_back(dst[4]);
	}

	sort(dists.begin(), dists.end());
	double mn = *min_element(dists.begin(), dists.end());
	double mx = *max_element(dists.begin(), dists.end());

	double threshold = -1.0;
	for (size_t j = 1; j < ncells; j++) {
		if (dists[j] - dists[j - 1] >= 2 * (mx - mn) / ((double) ncells)) {
			threshold = dists[j - 1];
			if (ind_thresh != NULL)
				*ind_thresh = j - 1;
			break;
		}
	}

	//Debug
	if (threshold == -1.0) {
		cerr << "Fail: couldn\'t find distance threshold\n";
	} else {
		cerr << "Distance threshold: " << threshold
				<< " corresponding to index " << *ind_thresh << "\n";
		cerr << 100.0 * (cells.size() - int(*ind_thresh)) / ((double) ncells)
				<< "\% of cells are noise\n";
	}

	return threshold;
}
