#include "cellset.h"

#include <iostream>
#include <ctime>
#include <algorithm>

using std::cerr;
using std::cin;
using std::cout;
using std::endl;
using std::min;

void RunAndPrintTime(void (CellSet::*f)(), CellSet *cs, string function_name) {
	time_t start, end;
	cerr << "Running " << function_name << "...\n";
	time(&start);
	(cs->*f)();
	time(&end);
	cerr << "Done. Elapsed time: " << end - start << "\n\n";
}

int main(int argc, char **argv) {
	cerr << "Program starded. Initializing CellSet object...\n";

	// Initializing CellSet object
	CellSet cs = CellSet(
	//Num neighbors for kNN graph
			min(cs.GetNumCells() - 1, (size_t) 20),

			//Exclude cells far away from all others
			true,

			//istream where data comes from
			cin,

			//ostream where data goes to
			cout
		);

	RunAndPrintTime(&CellSet::ReadFromStream, &cs, "read input from stdin");
	RunAndPrintTime(&CellSet::CreateVpTree, &cs, "vp tree");
	RunAndPrintTime(&CellSet::BuildNNGraph, &cs, "NN graph");
	RunAndPrintTime(&CellSet::CountConnectedComponents, &cs, "connected components before median");
	RunAndPrintTime(&CellSet::MedianFilter, &cs, "median filter");
	RunAndPrintTime(&CellSet::CountConnectedComponents, &cs, "count connected components after median");
	RunAndPrintTime(&CellSet::WriteToStream, &cs, "write to stdout");
		cerr << "Done!" << endl;
	return EXIT_SUCCESS;
}

