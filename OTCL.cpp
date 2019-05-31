#include <bits/stdc++.h>
#include "argparse.h"
using namespace std;
typedef vector<vector<int> > Graph;
typedef unordered_map<long long, long long> SharedEdge;
typedef unordered_map<long long, queue<int> > QueueSuperLeaf;
const int error = 1000;
const int size_default = 3;
int MAX_SIZE = error, MAX_LEAF = error;
int node = 0, s_node = 0, idx = 0, countTree = 0; // internal_node, subdividing_node
Graph pTree;
int dfsParser(string &newick, int parent) {
    int curNode = node;
    vector<int> adjacent;
    if(parent != 0) adjacent.push_back(parent);
    if(newick[idx] == '(') {
        curNode = node++;
        if(node == 1) node = countTree*MAX_LEAF+1;

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

/* ---------- find T|R and superleaf attach to edge ---------- */
void delete_node(Graph &T_R, int node, int del_node) {
    for(int i=0;i<T_R[node].size();i++) {
        if(T_R[node][i] == del_node) {
            T_R[node].erase(T_R[node].begin() + i);
            break;
        }
    }
}
bool DFS_FindT_R(int root, int parent, Graph &T, Graph &T_R, vector<int> &R, bool *isLeafR, Graph &superleaf) {
    if(isLeafR[root]) {
        return true;
    }
    for(auto child : T[root]) {
        if(child != parent) {
            if(DFS_FindT_R(child, root, T, T_R, R, isLeafR, superleaf) == true) {
                T_R[root].push_back(child);
                T_R[child].push_back(root);
            }
            else {
                superleaf[root].push_back(child);
                delete_node(T, child, root);
            }
        }
    }
    return T_R[root].empty() == false;
}
int DFS_Suppresed(int root, int parent, Graph &T_R, Graph &superleaf) {
    if(T_R[root].size() == 2) {
        for(auto S : superleaf[root]) {
            superleaf[parent].push_back(S);
        }

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
        return DFS_Suppresed(next_child, root, T_R, superleaf);
    }
    else {
        vector<int> child;
        for(auto c : T_R[root]) {
            if(c != parent) {
                child.push_back(c);
            }
        }
        for(auto c : child) {
            int childNode = DFS_Suppresed(c, root, T_R, superleaf);
            for(auto S : superleaf[root]) {
                superleaf[childNode].push_back(S);
            }
            superleaf[root].clear();
        }
        return root;
    }
}
void FindT_R(int r, Graph &T, Graph &T_R, vector<int> &R, vector<int> &Internal_T, bool *isLeafR, Graph &superleaf) {
    int root;
    for(auto node : Internal_T) {
        for(auto child : T[node]) {
            if(child == r)
                root = node;
        }
    }
    DFS_FindT_R(root, root, T, T_R, R, isLeafR, superleaf);

    root = r;
    DFS_Suppresed(root, root, T_R, superleaf);
}

/* ---------- Calculate Hash ----------*/
long long calculateEdgeNumber(int node1, int node2) {
    if(node1 < node2) return node1*1ll*(MAX_SIZE+1) + node2;
    return node2*1ll*(MAX_SIZE+1) + node1;
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

void DFS_SUPERLEAF(Graph &T, int u, int parent, long long e_hat, SharedEdge *sharedEdge, Graph &superleaf, QueueSuperLeaf *Q) {
    for(auto child : T[u]) {
        if(child != parent) {
            long long edgeNum = calculateEdgeNumber(u, child);
            long long e_tilde = e_hat;
            if(sharedEdge[1].find(edgeNum) != sharedEdge[1].end()) {
                e_tilde = edgeNum;

                // type I Superleaf
                for(auto S : superleaf[child]) {
                    Q[0][sharedEdge[1][e_tilde]].push(S);
                }
            }
            else {
                // type II Superleaf
                for(auto S : superleaf[child]) {
                    Q[1][sharedEdge[1][e_tilde]].push(S);
                }
            }
            DFS_SUPERLEAF(T, child, u, e_tilde, sharedEdge, superleaf, Q);
        }
    }
}


void DFS_AddSubTree(Graph &t, Graph &T, int u, int parent, int curNode_t, bool *isLeafS) {
    for(auto child : T[u]) {
        if(child != parent){
            if(isLeafS[child] == true) {
                t[curNode_t].push_back(child);
                t[child].push_back(curNode_t);
            }
            else {
                t[curNode_t].push_back(s_node);
                t[s_node].push_back(curNode_t);

                s_node++;
                DFS_AddSubTree(t, T, child, u, s_node-1, isLeafS);
            }
        }
    }
}
void DFS_AttachSuperleaf(Graph &t, Graph &T, int u, int parent, QueueSuperLeaf *Q, bool *isLeafS) {
    vector<int> child;
    for(auto c : t[u]) {
        if(c != parent) {
            child.push_back(c);
        }   
    }
    for(auto c : child) {
        DFS_AttachSuperleaf(t, T, c, u, Q, isLeafS);
    }

    int v = parent;
    long long edgeNum = calculateEdgeNumber(v, u);
    if(Q[0].find(edgeNum) != Q[0].end()) {
        while(Q[0][edgeNum].empty() == false) {
            int superleaf = Q[0][edgeNum].front();
            Q[0][edgeNum].pop();

            delete_node(t, v, u);
            delete_node(t, u, v);
            t[u].push_back(s_node);
            t[v].push_back(s_node);
            t[superleaf].push_back(s_node);

            t[s_node].push_back(superleaf);
            t[s_node].push_back(u);
            t[s_node].push_back(v);

            v = s_node;
            s_node++;

            DFS_AddSubTree(t, T, superleaf, s_node, superleaf, isLeafS);
        }
    }

    if(Q[1].find(edgeNum) != Q[1].end()) {
        while(Q[1][edgeNum].empty() == false) {
            int superleaf = Q[1][edgeNum].front();
            Q[1][edgeNum].pop();

            delete_node(t, v, u);
            delete_node(t, u, v);
            t[u].push_back(s_node);
            t[v].push_back(s_node);
            t[superleaf].push_back(s_node);

            t[s_node].push_back(superleaf);
            t[s_node].push_back(u);
            t[s_node].push_back(v);

            u = s_node;
            s_node++;

            DFS_AddSubTree(t, T, superleaf, s_node, superleaf, isLeafS);
        }
    }
}

/* ---------- Algorithm ---------- */   
Graph Algorithm(Graph &T, Graph &t) {
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
        int r = 1;

        // find T|R
        Graph superleaf;
        superleaf.resize(MAX_SIZE);

        Graph T_R;
        T_R.resize(MAX_SIZE);
        FindT_R(r, T, T_R, R, Internal_T, isLeafR, superleaf);

        // calculate shared edge
        long long hashNode[MAX_LEAF+10];
        for(auto leaf : S) {
            hashNode[leaf] = (rand()*(1ll<<30)) + rand()%(1ll<<30);
        }
        
        SharedEdge sharedEdge[2];
        calculate_Bipartition(T_R, t, hashNode, R, isLeafR, sharedEdge);

        long long e_hat = -1;
        QueueSuperLeaf Q[2];

        DFS_SUPERLEAF(T_R, r, r, e_hat, sharedEdge, superleaf, Q);

        // attach each superleaf
        Graph _t = t;
        DFS_AttachSuperleaf(_t, T, r, r, Q, isLeafS);

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
    parser.add_argument("--seed", "an integer", false);
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

    if(seed == 0) srand(time(0));
    else srand(seed);
    if(m == 0) m = size_default + error;
    else m += error;
    if(n == 0) n = size_default + error;
    else n += error;

    MAX_SIZE = 4*n;
    MAX_LEAF = n;
    s_node = 3*n+1;

    string newick_T, newick_t;
    ifstream file;
    file.open(file_T);
    while(file) getline(file, newick_T);
    file.close();

    file.open(file_t);
    while(file) getline(file, newick_t);
    file.close();

    Graph T, t;
    idx = 0, node = 0, countTree = 1;
    T = Parser(newick_T);

    idx = 0, node = 0, countTree = 2;
    t = Parser(newick_t);

    Graph result = Algorithm(T, t);
    int root;
    for(int i=MAX_LEAF+1;i<MAX_SIZE;i++) {
        if(result[i].empty() == false) {
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