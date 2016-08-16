#include "cellset.h"

#include <iostream>
#include <ctime>
#include <algorithm>

using std::cerr;
using std::cin;
using std::endl;
using std::min;

//Remove to disable verbose
#define DEBUG true


void RunAndPrintTime(void (CellSet::*f)(), CellSet *cs, string function_name) {
	time_t start,end;
	cerr << "Running " << function_name << "...\n";
	time(&start);
	(cs->*f)();
	time(&end);
	cerr << "Finished " << function_name << "\n";
	cerr << "Elapsed time: " << end-start << "\n\n";
}

int main(int argc, char **argv) {
	cerr << "Program starded. Initializing CellSet object...\n";

	// Initializing CellSet object
	CellSet cs = CellSet(
			//Num neighbors for kNN graph
			min(cs.GetNumCells(), (size_t)20),

			//Exclude cells far away from all others
			true,

			//istream where data comes from
			cin
		);

	RunAndPrintTime(&CellSet::ReadFromStream, &cs, "read input from stdin");
	RunAndPrintTime(&CellSet::CreateVpTree, &cs, "create vp tree");
	RunAndPrintTime(&CellSet::BuildNNGraph, &cs, "build NN graph");
	RunAndPrintTime(&CellSet::EstimateReads, &cs, "apply median filter");

	cerr << "Done!" << endl;
	return EXIT_SUCCESS;
}

