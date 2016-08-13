#include "cell.h"
#include <algorithm>
using std::cerr;
using std::min;

// Class that will store all useful data about the cell:
// nonzero coordinates, cluster index, name and index
Cell::Cell(string _name, size_t _index) : name(_name), index(_index) {}

void Cell::insert(int indv, double valv) {
	ind.push_back(indv);
	val.push_back(valv);
}

void Cell::print() {
	int sz = ind.size();
	cerr << name << "\t";
	cerr << sz << "\t";
	for (int i = 0; i < min(sz, 10); i++) {
		cerr << ind[i] << " " << val[i] << "\t";
	}
	cerr << "\n";
}
