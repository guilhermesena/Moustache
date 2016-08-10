#include <iostream>
#include <cstdio>
#include <vector>
#include <sstream>
#include <string>
#include <algorithm>
#include <ctime>
#define eps 1e-6

using std::string;
using std::vector;
using std::getline;
using std::cin;
using std::cout;
using std::istringstream;
using std::getline;
using std::ios;
using std::min;

class cell_data {
public:
    vector<size_t> ind;
    vector<double> val;
    string name;

    cell_data(string _name) {
        name = _name;
    }

    void insert(int indv, double valv) {
        ind.push_back(indv);
        val.push_back(valv);
    }

    void print(){
        int sz = ind.size();
        cout << name << "\t";
        cout << sz << "\t";
        for(int i = 0; i < min(sz,10); i++) {
            cout << ind[i] << " " << val[i] << "\t";
        }
        cout << "\n";
    }
};   

static double dist(cell_data &ca, cell_data &cb) {
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

    return ans;
}


int main(int argc, char **argv) {
    string cellname, line;
    char genename[100];
    float count;
    size_t ncells = 0, gene = 0, i;
    time_t start, read_input, pairwise_distances;
    time(&start);

    cout << "Counting cells...\n";
    getline(cin,line);
    istringstream rownames(line);
    vector<cell_data> cells;
    while(rownames >> cellname) {
        cells.push_back(cell_data(cellname));
    }

    ncells = cells.size();

    cout << "Reading gene data...\n";
    while (scanf(" %s", genename) > 0) {
        if((gene%1000) == 0) cout << "Reading gene " << genename << "\n";
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
    cout << "Time to read input: " << read_input - start << "\n";

    cout << "total cells: " << cells.size() << "\n";
    cout << "total genes: " << gene << "\n";
    for (int i = 0; i < min((size_t)20, ncells); i++) cells[i].print();

    cout << "Computing pairwise distances...\n";
    for (int i = 0; i < ncells; i++) {
        for(int j = 0; j < i; j++) {
            double d = dist(cells[i], cells[j]);
            cout << "Distance between " << i << "and " << j << "is " << d << "\n";
        }
    }
    time (&pairwise_distances);
    cout << "Time to find pwd: " << pairwise_distances - read_input << "\n";
    cout << "Done!\n";
}
