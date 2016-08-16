#include "cell.h"
#include "vptree.h"

#include <utility>
#include <istream>
#include <algorithm>
#include <unordered_map>

//Debug
#include <iostream>

using std::pair;
using std::unordered_map;
using std::istream;
using std::min;
using std::max;
using std::make_pair;

//Hashing of pairs to use in unordered_map
template<class T>
inline void hash_combine(std::size_t & seed, const T & v) {
	std::hash<T> hasher;
	seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

namespace std {
template<typename S, typename T> struct hash<pair<S, T>> {
	inline size_t operator()(const pair<S, T> & v) const {
		size_t seed = 0;
		::hash_combine(seed, v.first);
		::hash_combine(seed, v.second);
		return seed;
	}
};
}

//A class that will read and store the cell data
class CellSet {
public:
	CellSet();
	~CellSet();
	void ReadFromStream(istream&);
	void CreateVpTree();
	void BuildNNGraph(int NumNeighbors, bool ExcludeNoise, size_t *ind_thresh =
			NULL);
	void EstimateReads();

	static unordered_map<pair<size_t, size_t>, double> dist_memory;

private:
	struct NNPair {
		size_t CellIndex;
		double distance;

		NNPair(size_t _CellIndex, double _distance) :
				CellIndex(_CellIndex), distance(_distance) {
		}

		bool operator <(const NNPair &rhs) const {
			return distance < rhs.distance;
		}
	};

	vector<Cell*> cells;

	// Static function to calculate euclidean distances between cells
	static double dist(Cell* ca, Cell* cb) {
		size_t mn = min(ca->index, cb->index);
		size_t mx = max(ca->index, cb->index);
		if (dist_memory[make_pair(mn, mx)]) {
			return dist_memory[make_pair(mn, mx)];
		}
		double ans = 0;
		int pa, pb;
		size_t sa = ca->ind.size(), sb = cb->ind.size();

		for (pa = 0, pb = 0; pa < sa && pb < sb;) {
			if (ca->ind[pa] == cb->ind[pb]) {
				ans += (ca->val[pa] - cb->val[pb])
						* (ca->val[pa] - cb->val[pb]);
				pa++;
				pb++;
			} else if (ca->ind[pa] < cb->ind[pb]) {
				ans += ca->val[pa] * ca->val[pa];
				pa++;
			} else {
				ans += cb->val[pb] * cb->val[pb];
				pb++;
			}
		}
		while (pa < sa) {
			ans += ca->val[pa] * ca->val[pa];
			pa++;
		}
		while (pb < sb) {
			ans += cb->val[pb] * cb->val[pb];
			pb++;
		}

		return dist_memory[make_pair(mn, mx)] = sqrt(ans);
	}
	VpTree<Cell*, dist> vpt;
};
