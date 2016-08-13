#include "cell.h"
#include "vptree.h"

#include <utility>
#include <istream>
#include <algorithm>
#include <unordered_map>
#include <iostream>

//Debug
#include <iostream>


using std::pair;
using std::unordered_map;
using std::istream;
using std::min;
using std::max;
using std::make_pair;
using std::cerr;

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
	double GetDistanceThreshold(size_t *ind_thresh = NULL);
	static unordered_map<pair<size_t, size_t>, double> dist_memory;

private:
	vector<Cell> cells;
	// Static function to calculate euclidean distances between cells
	static double dist(const Cell &ca, const Cell &cb) {
		size_t mn = min(ca.index, cb.index);
		size_t mx = max(ca.index, cb.index);
		if (dist_memory[make_pair(mn, mx)]) {
			std::cerr << "distance calculation avoided!\n";
			return dist_memory[make_pair(mn, mx)];
		}
		double ans = 0;
		int pa = 0, pb = 0;
		int sa = ca.ind.size(), sb = cb.ind.size();

		for (; pa < sa && pb < sb;) {
			if (ca.ind[pa] == cb.ind[pb]) {
				ans += (ca.val[pa] - cb.val[pb]) * (ca.val[pa] - cb.val[pb]);
				pa++;
				pb++;
			} else if (ca.ind[pa] < cb.ind[pb]) {
				ans += ca.val[pa] * ca.val[pa];
				pa++;
			} else {
				ans += cb.val[pb] * cb.val[pb];
				pb++;
			}
		}
		while (pa < sa) {
			ans += ca.val[pa] * ca.val[pa];
			pa++;
		}
		while (pb < sb) {
			ans += cb.val[pb] * cb.val[pb];
			pb++;
		}

		return dist_memory[make_pair(mn, mx)] = sqrt(ans);
	}
	VpTree<Cell, dist> vpt;
};
