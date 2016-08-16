#include "cellset.h"

#include <iostream>
#include <ctime>

using std::cerr;
using std::cin;

//Remove to disable verbose
#define DEBUG true

int main(int argc, char **argv) {
	time_t time_start, time_read_input, time_vp_tree, time_dists;
	cerr << "Program starded. Reading input...\n";
	time(&time_start);

	// Initializing CellSet object
	CellSet cs = CellSet();

	// Read dataset from stdin (change cin to ifstream if necessary)
	cs.ReadFromStream(cin);
	time(&time_read_input);
	cerr << "Time to read input: " << time_read_input - time_start
			<< " seconds\n";

	//VP tree construction
	cerr << "Creating VPtree...\n";
	cs.CreateVpTree();
	time(&time_vp_tree);
	cerr << "Time to build vp tree: " << time_vp_tree - time_read_input
			<< " seconds\n";

	//Create KNN graph
	cerr << "Creating kNN graph...\n";
	size_t ind_threshold;
	cs.BuildNNGraph(20, true, &ind_threshold);
	time(&time_dists);
	cerr << "Time to compute distances: " << time_dists - time_vp_tree
			<< " seconds \n";

	//Denoise
	cerr << "Estimating true count for reads\n";
	cs.EstimateReads();

	cerr << "Done!" << endl;
	return EXIT_SUCCESS;
}

