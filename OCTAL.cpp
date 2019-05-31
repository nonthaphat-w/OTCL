#include <bits/stdc++.h>
#include "argparse.h"
using namespace std;
typedef vector<vector<int> > Graph;
typedef unordered_map<long long, long long> SharedEdge;
const int error = 1000;
const int size_default = 3;
int MAX_SIZE = error, MAX_LEAF = error;
int node = 0, s_node = 0, idx = 0; // internal_node, subdividing_node
Graph pTree;
int dfsParser(string &newick, int parent) {
    int curNode = node;
    vector<int> adjacent;
    if(parent != 0) adjacent.push_back(parent);
    if(newick[idx] == '(') {
        curNode = node++;
        if(node == 1) node = MAX_LEAF+1;

        // left child (
        idx++;
        adjacent.push_back(dfsParser(newick, curNode));

        // right child ,
        idx++;
        adjacent.push_back(dfsParser(newick, curNode));

        // )
        idx++;
    }
    else {
        int leaf_node = 0;
        while(newick[idx] != ',' && newick[idx] != ')') {
            leaf_node = leaf_node*10 + newick[idx]-'0';
            idx++;
        }
        curNode = leaf_node;
    }
    pTree[curNode] = adjacent;
    if(curNode == 0) {
        pTree[pTree[0][0]].push_back(pTree[0][1]);
        pTree[pTree[0][1]].push_back(pTree[0][0]);
    }
    pTree[0] = vector<int>();
    return curNode;
}
Graph Parser(string newick) {
    // remove space
    newick.erase(remove_if(newick.begin(), newick.end(), ::isspace), newick.end());

    pTree.clear();
    pTree.resize(MAX_SIZE);
    idx = 0, node = 0;

    // parse tree
    dfsParser(newick, 0);
    return pTree;
}

/* ---------- find T|R ---------- */
bool DFS_FindT_R(int root, int parent, Graph &T, Graph &T_R, vector<int> &R, bool *isLeafR) {
    if(isLeafR[root]) {
        return true;
    }
    for(auto child : T[root]) {
        if(child != parent) {
            if(DFS_FindT_R(child, root, T, T_R, R, isLeafR) == true) {
                T_R[root].push_back(child);
                T_R[child].push_back(root);
            }
        }
    }
    return !T_R[root].empty();
}
void delete_node(Graph &T_R, int node, int del_node) {
    for(int i=0;i<T_R[node].size();i++) {
        if(T_R[node][i] == del_node) {
            T_R[node].erase(T_R[node].begin() + i);
            break;
        }
    }
}
void DFS_Suppresed(int root, int parent, Graph &T_R) {
    if(T_R[root].size() == 2) {
        int next_child;
        for(auto child : T_R[root]) {
            if(child != parent) {
                next_child = child;
            }
        }
        T_R[root] = vector<int>();
        delete_node(T_R, parent, root);
        delete_node(T_R, next_child, root);

        root = parent;
        T_R[root].push_back(next_child);
        T_R[next_child].push_back(root);
        DFS_Suppresed(next_child, root, T_R);
    }
    else {
        vector<int> child;
        for(auto c : T_R[root]) {
            if(c != parent) {
                child.push_back(c);
            }
        }
        for(auto c : child) {
            DFS_Suppresed(c, root, T_R);
        }
    }
}
void FindT_R(Graph &T, Graph &T_R, vector<int> &R, vector<int> &Internal_T, bool *isLeafR) {
    int root;
    for(auto node : Internal_T) {
        for(auto child : T[node]) {
            if(isLeafR[child] == true)
                root = node;
        }
    }
    DFS_FindT_R(root, root, T, T_R, R, isLeafR);

    root = R[0];
    DFS_Suppresed(root, root, T_R);
}

/* ---------- Calculate Hash ----------*/
long long calculateEdgeNumber(int node1, int node2) {
    return node1*1ll*(MAX_SIZE+1) + node2;
}
long long DFS_calHash(int root, int parent, Graph &T, long long *hashNode, unordered_map<long long, long long> &Bipartition_T, bool *isLeafR) {
    long long hashValue = 0;
    for(auto child : T[root]) {
        if(child != parent) {
            hashValue ^= DFS_calHash(child, root, T, hashNode, Bipartition_T, isLeafR);
        }
    }

    if(isLeafR[root] == true) {
        hashValue ^= hashNode[root];
    }
    if(root != parent) {
        long long edgeNum = 0;
        if(parent < root) edgeNum = calculateEdgeNumber(parent, root);
        else edgeNum = calculateEdgeNumber(root, parent);
        Bipartition_T[edgeNum] = hashValue;
    }

    return hashValue;
}
void calculate_Bipartition(Graph &T, Graph &t, long long *hashNode, vector<int> &R, bool *isLeafR, SharedEdge *sharedEdge) {
    sharedEdge[0].clear();
    sharedEdge[1].clear();
    long long totalHash = 0;
    for(auto leaf : R) totalHash ^= hashNode[leaf];

    unordered_map<long long, long long> Bipartition_T, Bipartition_t;
    int root = R[0];
    DFS_calHash(root, root, T, hashNode, Bipartition_T, isLeafR);
    DFS_calHash(root, root, t, hashNode, Bipartition_t, isLeafR);

    SharedEdge edgeHash;
    for(auto edge : Bipartition_T) {
        edgeHash[edge.second] = edge.first;
    }
    for(auto edge : Bipartition_t) {
        if(edgeHash.find(edge.second^totalHash) != edgeHash.end()) edge.second ^= totalHash;
        if(edgeHash.find(edge.second) != edgeHash.end()) {
            sharedEdge[0][edge.first] = edgeHash[edge.second];
            sharedEdge[1][edgeHash[edge.second]] = edge.first;
        }
    }
}

/* ---------- ADD LEAF ----------*/
bool DFS_Find_y(int root, int parent, int y, Graph &T, vector<int> &path) {
    if(root == y) {
        path.push_back(y);
        return true;
    }
    else {
        for(auto child : T[root]) {
            if(child != parent) {
                if(DFS_Find_y(child, root, y, T, path) == true) {
                    path.push_back(root);
                    return true;
                }
            }
        }
    }
    return false;
}

void ADDLEAF(int x, int y, Graph &t, Graph &T, SharedEdge *sharedEdge) {
    // root at v
    int v = T[x][0];

    // find path of v to y
    vector<int> path;
    DFS_Find_y(v, x, y, T, path);
    reverse(path.begin(), path.end());
    bool ok = 0;
    for(auto child : T[v]) {
        if(child != x && child != path[1]) {
            path[0] = child;
            break;
        }
    }

    // find e' in t defining the same biparition as e
    int previous_node = -1;
    for(auto node : path) {
        if(previous_node != -1) {
            long long edgeNum_T = 0;
            if(previous_node < node) edgeNum_T = calculateEdgeNumber(previous_node, node);
            else edgeNum_T = calculateEdgeNumber(node, previous_node);

            if(sharedEdge[1].find(edgeNum_T) != sharedEdge[1].end()) {
                long long edgeNum_t = sharedEdge[1][edgeNum_T];

                // e' = (node1, node2)
                int node1 = edgeNum_t/(MAX_SIZE+1), node2 = edgeNum_t%(MAX_SIZE+1);
                delete_node(t, node1, node2);
                delete_node(t, node2, node1);

                // attach x to e'
                t[node1].push_back(s_node);
                t[node2].push_back(s_node);
                t[x].push_back(s_node);

                t[s_node].push_back(node1);
                t[s_node].push_back(node2);
                t[s_node].push_back(x);

                s_node++;
                return;
            }
        }
        previous_node = node;
    }
}

Graph OCTAL(Graph &T, Graph &t) {
    bool isLeafR[MAX_SIZE+10], isLeafS[MAX_SIZE+10];
    memset(isLeafS, 0, sizeof(isLeafS));
    memset(isLeafR, 0, sizeof(isLeafR));
    vector<int> R, S;
    vector<int> Internal_T, Internal_t;
    for(int i=1;i<=MAX_LEAF;i++){
        if(T[i].size() == 1) {
            S.push_back(i);
            isLeafS[i] = true;
        }
        if(t[i].size() == 1) {
            R.push_back(i);
            isLeafR[i] = true;
        }
    }
    for(int i=MAX_LEAF+1;i<MAX_SIZE;i++){
        if(T[i].empty() == false) {
            Internal_T.push_back(i);
        }
        if(t[i].empty() == false) {
            Internal_t.push_back(i);
        }
    }

    if(S == R) {
        return t;
    }
    else {
        // find T|R
        Graph T_R;
        T_R.resize(MAX_SIZE);
        FindT_R(T, T_R, R, Internal_T, isLeafR);

        // calculate shared edge
        long long hashNode[MAX_LEAF+10];
        for(auto leaf : S) {
            hashNode[leaf] = (rand()*(1ll<<30)) + (rand()%(1ll<<30));
        }
        
        SharedEdge sharedEdge[2];
        calculate_Bipartition(T_R, t, hashNode, R, isLeafR, sharedEdge);

        vector<int> _R = R;
        Graph _t = t;
        Graph _T;
        _T.resize(MAX_SIZE);

        // S\R
        for(int x=R.size()+1;x<=S.size();x++){
            _R.push_back(x);
            isLeafR[x] = true;

            for(auto node : _R)
                _T[node].clear();
            for(auto node : Internal_T)
                _T[node].clear();
            FindT_R(T, _T, _R, Internal_T, isLeafR);

            int y = 1;
            ADDLEAF(x, y, _t, _T, sharedEdge);
            calculate_Bipartition(_T, _t, hashNode, _R, isLeafR, sharedEdge);
        }
        return _t;
    }
}

string encodeNewick(int node, int parent, int root, Graph &T) {
    string result = "";

    if(T[node].size() == 1) {
        result += to_string(node);
    }
    else if(node == root) {
        int count_child = 0;
        result += "((";
        for(auto child : T[node]) {
            if(child != parent) {
                result += encodeNewick(child, node, root, T);
                if(count_child == 0) result += ", ";
                else if(count_child == 1) result += "), ";
                else result += ")";
                count_child++;
            }
        }
        result += ";";
    }
    else {
        int count_child = 0;
        result += "(";
        for(auto child : T[node]) {
            if(child != parent) {
                result += encodeNewick(child, node, root, T);
                if(count_child == 0) result += ", ";
                else result += ")";
                count_child++;
            }
        }
    }
    return result;
}

// input_T, input_t, seed, size of T, size of t
int main(int argc, char *argv[]){
    ArgumentParser parser("Argument parser example");
    parser.add_argument("-i", "input files", true);
    parser.add_argument("-b", "input files", true);
    parser.add_argument("-seed", "an integer", false);
    parser.add_argument("-m", "an integer", false);
    parser.add_argument("-n", "an integer", false);
    parser.add_argument("-o", "input files", true);

    try {
        parser.parse(argc, argv);
    } catch (const ArgumentParser::ArgumentNotFound& ex) {
        std::cout << ex.what() << std::endl;
        return 0;
    }
    if (parser.is_help()) return 0;

    auto file_t = parser.get<string>("i");
    auto file_T = parser.get<string>("b");
    int seed = parser.get<int>("seed");
    int m = parser.get<int>("m");
    int n = parser.get<int>("n");
    auto file_O = parser.get<string>("o");
    if(seed == 1) srand(time(0));
    else srand(seed);
    if(m == 1) m = size_default + error;
    else m += error;
    if(n == 1) n = size_default + error;
    else n += error;

    MAX_SIZE = 3*n;
    MAX_LEAF = n;
    s_node = 2*n+1;

    string newick_T, newick_t;
    ifstream file;
    file.open(file_T);
    while(file) getline(file, newick_T);
    file.close();

    file.open(file_t);
    while(file) getline(file, newick_t);
    file.close();

    Graph T, t;
    idx = 0;
    T = Parser(newick_T);
    idx = 0;
    t = Parser(newick_t);

    Graph result = OCTAL(T, t);
    int root;
    for(int i=MAX_LEAF+1;i<MAX_SIZE;i++) {
        if(!result[i].empty()) {
            root = i;
            break;
        }
    }

    string newick_result = encodeNewick(root, root, root, result);
    ofstream write;
    write.open(file_O);
    write << newick_result;
    write.close();
}