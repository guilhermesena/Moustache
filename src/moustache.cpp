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
using std::priority_queue;
using std::reverse;

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

class cell {
public:
    vector<size_t> ind;
    vector<double> val;
    string name;

    cell(string _name) {
        name = _name;
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

double dist(const cell &ca, const cell &cb) {
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
}


int main(int argc, char **argv) {
    string cellname, line;
    char genename[100];
    float count;
    size_t ncells = 0, gene = 0, i;
    time_t start, read_input, vp_tree, dists;
    time(&start);

    cerr << "Counting cells...\n";
    getline(cin,line);
    istringstream rownames(line);
    vector<cell> cells;
    while(rownames >> cellname) {
        cells.push_back(cell(cellname));
    }

    ncells = cells.size();

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

    time(&read_input);
    cerr << "Time to read input: " << read_input - start << " seconds\n";

    cerr << "total cells: " << cells.size() << "\n";
    cerr << "total genes: " << gene << "\n";
    for (int i = 0; i < min((size_t)20, ncells); i++) cells[i].print();

    cerr << "Making VPtree...\n";
    VpTree<cell, dist> vpt = VpTree<cell, dist>();
    vpt.create(&cells);
    time(&vp_tree);
    cerr << "Time to build vp tree: " << vp_tree - read_input << " seconds\n";
    for (int i = 0; i < min((size_t)20, ncells); i++) cells[i].print();

    cerr << "Computing 4th distances...\n";
    for (int i = 0; i < min(ncells,(size_t)500); i++) {
        vector<cell> res;
        vector<double> dst;
        vpt.search(cells[i], 6, &res, &dst);
        cout << dst[4] << " ";
    }
    time(&dists);
    cerr << "Time to compute distances: " << dists - vp_tree << " seconds \n";
}

