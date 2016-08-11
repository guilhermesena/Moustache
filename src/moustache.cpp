#include <iostream>
#include <cstdio>
#include <vector>
#include <sstream>
#include <string>
#include <algorithm>
#include <ctime>
#include <cmath>
#include <queue>
#include <limits> 
#include <cstdlib>
#include <unordered_map>
#include <utility>
#define eps 1e-6

using std::queue;
using std::string;
using std::vector;
using std::getline;
using std::cin;
using std::cout;
using std::cerr;
using std::istringstream;
using std::getline;
using std::ios;
using std::min;
using std::max;
using std::priority_queue;
using std::reverse;
using std::sort;
using std::min_element;
using std::max_element;
using std::unordered_map;
using std::pair;
using std::make_pair;

//map to avoid recalculating distance between cells
template <class T>
inline void hash_combine(std::size_t & seed, const T & v) {
      std::hash<T> hasher;
        seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

namespace std{
    template<typename S, typename T> struct hash<pair<S, T>> {
        inline size_t operator()(const pair<S, T> & v) const {
            size_t seed = 0;
            ::hash_combine(seed, v.first);
            ::hash_combine(seed, v.second);
            return seed;
        }
    };
}

unordered_map<pair<size_t, size_t>, double> dist_memory;

// Class that will store all useful data about the cell:
// nonzero coordinates, cluster index, name and index
class cell {
public:
    vector<size_t> ind;
    vector<double> val;
    string name;
    int index;


    cell(string _name, size_t _index) {
        name = _name;
        index = _index;
    }

    void insert(int indv, double valv) {
        ind.push_back(indv);
        val.push_back(valv);
    }

    void print(){
        int sz = ind.size();
        cerr << name << "\t";
        cerr << sz << "\t";
        for(int i = 0; i < min(sz,10); i++) {
            cerr << ind[i] << " " << val[i] << "\t";
        }
        cerr << "\n";
    }
};   

// Static function to calculate euclidean distances between cells
static double dist(const cell &ca, const cell &cb) {
    size_t mn = min(ca.index, cb.index);
    size_t mx = max(ca.index, cb.index);
    if(dist_memory[make_pair(mn,mx)]){
        cerr << "Redundant computation avoided!\n";
        return dist_memory[make_pair(mn,mx)];
    }
    double ans = 0;
    int pa = 0, pb = 0;
    int sa = ca.ind.size(), sb = cb.ind.size();

    for(; pa < sa && pb < sb;) {
        if(ca.ind[pa] == cb.ind[pb]) {
            ans += (ca.val[pa] - cb.val[pb])*(ca.val[pa]-cb.val[pb]);
            pa++;
            pb++;
        } else if(ca.ind[pa] < cb.ind[pb]) {
            ans += ca.val[pa]*ca.val[pa];
            pa++;
        } else {
            ans += cb.val[pb]*cb.val[pb];
            pb++;
        }
    }
    while(pa < sa) {
        ans += ca.val[pa]*ca.val[pa];
        pa++;
    }
    while(pb < sb) {
        ans += cb.val[pb]*cb.val[pb];
        pb++;
    }

    return sqrt(ans);
    return dist_memory[make_pair(mn,mx)] = sqrt(ans);
}


// VP tree data structure for fast retrieval of nearest neighbors
template<typename T, double (*distance)( const T&, const T& )>
class VpTree {
public:
    VpTree() : _root(0) {}
    ~VpTree() {
        delete _root;
    }
    
    void create( vector<T>* items ) {
        delete _root;
        _items = items;
        _root = buildFromPoints(0, items->size());
    }
    
    void search( const T& target, int k, vector<T>* results, vector<double>* distances) {
        priority_queue<HeapItem> heap;
        _tau = std::numeric_limits<double>::max();
        search( _root, target, k, heap );
        results->clear(); distances->clear();
        while( !heap.empty() ) {
            results->push_back( (*_items)[heap.top().index] );
            distances->push_back( heap.top().dist );
            heap.pop();
        }
        
        reverse( results->begin(), results->end() );
        reverse( distances->begin(), distances->end() );
    }

private:
    vector<T> *_items;
    double _tau;
    struct Node {
        int index;
        double threshold;
        Node* left;
        Node* right;
        Node() :index(0), threshold(0.), left(0), right(0) {}
        
        ~Node() {
            delete left;
            delete right;
        }
    } * _root;
    
    struct HeapItem {
        HeapItem( int index, double dist) : index(index), dist(dist) {}
        
        int index;
        double dist;
        bool operator<( const HeapItem& o ) const {
            return dist < o.dist;   
        }
    };
    
    struct DistanceComparator {
        const T& item;
        DistanceComparator( const T& item ) : item(item) {}
        bool operator()(const T& a, const T& b) {
            return distance( item, a ) < distance( item, b );
        }
    };
    
    Node* buildFromPoints( int lower, int upper ) {
        if ( upper == lower ) {
            return NULL;
        }
        
        Node* node = new Node();
        node->index = lower;
        if ( upper - lower > 1 ) {
            int i = (int)((double)rand() / RAND_MAX * (upper - lower - 1) ) + lower;
            std::swap( (*_items)[lower], (*_items)[i] );
            int median = ( upper + lower ) / 2;
            
            std::nth_element(
                    _items->begin() + lower + 1, 
                    _items->begin() + median,
                    _items->begin() + upper,
                    DistanceComparator( (*_items)[lower] ));
            
            node->threshold = distance( (*_items)[lower], (*_items)[median] );
            node->index = lower;
            node->left = buildFromPoints( lower + 1, median );
            node->right = buildFromPoints( median, upper );
        }
        return node;
    }
    
    void search( Node* node, const T& target, int k, std::priority_queue<HeapItem>& heap ) {
        if ( node == NULL ) return;
        double dist = distance( (*_items)[node->index], target );
        if ( dist < _tau ) {
            if ( heap.size() == k ) heap.pop();
            heap.push( HeapItem(node->index, dist) );
            
            if ( heap.size() == k ) _tau = heap.top().dist;
        }
        
        if ( node->left == NULL && node->right == NULL ) {
            return;
        }
        
        if ( dist < node->threshold ) {
            if ( dist - _tau <= node->threshold ) {
                search( node->left, target, k, heap );
            }
            
            if ( dist + _tau >= node->threshold ) {
                search( node->right, target, k, heap );
            }
        
        } else {
            if ( dist + _tau >= node->threshold ) {
                search( node->right, target, k, heap );
            }
            
            if ( dist - _tau <= node->threshold ) {
                search( node->left, target, k, heap );
            }
        }
    }
};

static double get_distance_threshold (vector<cell> &cells, VpTree<cell, dist> &vpt, size_t *ind_thresh = NULL) {
    vector<double> dists;
    size_t ncells = cells.size();

    for (int i = 0; i < ncells; i++) {
        vector<cell> res;
        vector<double> dst;
        vpt.search(cells[i], 5, &res, &dst);

        //Debug
        cerr << cells[i].name << " ";
        for (size_t d = 0; d < dst.size(); d++) {
            cerr << dst[d] << " ";
        }
        cerr <<"\n";

        dists.push_back(dst[4]);
    }

    sort(dists.begin(), dists.end());
    double mn = *min_element(dists.begin(), dists.end());
    double mx = *max_element(dists.begin(), dists.end());

    double threshold = -1.0;
    size_t ind_threshold = -1;
    for(size_t j = 1; j < ncells; j++) {
        if(dists[j] - dists[j-1] >= 2*(mx - mn)/((double)ncells)){
            threshold = dists[j-1];
            if(ind_thresh != NULL) 
                *ind_thresh = j-1;
            break;
        }
    }
    return threshold;
}

 

int main(int argc, char **argv) {
    string cellname, line;
    char genename[100];
    float count;
    size_t ncells = 0, gene = 0, i;
    time_t time_start, time_read_input, time_vp_tree, time_dists;
    time(&time_start);

    cerr << "Counting cells...\n";
    getline(cin,line);
    istringstream rownames(line);
    vector<cell> cells;
    while(rownames >> cellname) {
        cells.push_back(cell(cellname, (int) ncells++));
    }


    cerr << "Reading gene data...\n";
    while (scanf(" %s", genename) > 0) {
        i = 0;
        while(i < ncells) {
            scanf(" %f", &count);
            if(count > eps)
                cells[i].insert(gene, count);
            i++;
        }
        gene++;
    }    

    time(&time_read_input);
    cerr << "Time to read input: " << time_read_input - time_start << " seconds\n";

    cerr << "total cells: " << cells.size() << "\n";
    cerr << "total genes: " << gene << "\n";
    for (int i = 0; i < min((size_t)20, ncells); i++) cells[i].print();

    cerr << "Making VPtree...\n";
    VpTree<cell, dist> vpt = VpTree<cell, dist>();
    vpt.create(&cells);
    time(&time_vp_tree);
    cerr << "Time to build vp tree: " << time_vp_tree - time_read_input << " seconds\n";
    for (int i = 0; i < min((size_t)20, ncells); i++) cells[i].print();

    cerr << "Computing 4th distances...\n";

    size_t ind_threshold;
    double threshold = get_distance_threshold(cells, vpt, &ind_threshold);

    if(threshold == -1.0) {
        cerr << "Fail: couldnt find distance threshold\n";
        return EXIT_FAILURE;
    } else {
        cerr << "Distance threshold: " << threshold << " corresponding to index " << ind_threshold << "\n";
        cerr << 100.0*(ncells - ind_threshold)/((double)ncells) << "\% of cells are noise\n";
    }

    time(&time_dists);
    cerr << "Time to compute distances: " << time_dists - time_vp_tree << " seconds \n";
    return EXIT_SUCCESS;
}


