//Zeynep Ozgur Gun, 29502

#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <chrono>

using namespace std;

const int N = 2;
static vector<int> rts;
// Stores the vertices of the original graph
int store[N + 1];

// Graph
vector<vector<int>> graph(N + 1, vector<int>(N + 1));

// Degree of the vertices
int d[N + 1];

// Max clique sizes for each vertex (for the heuristic algorithm)
int maxC[N + 1];


//Stores the maximum clique vertices for each clique

vector<int> BFclique;

vector<int> Hclique;


int randomInt(int min, int max) {
    return min + rand() % (max - min + 1);
}

void generateRandomGraph(int numNodes) {

    vector<int> all;
    for (int i = 0; i < numNodes; i++) {
        all.push_back(i + 1);
    }

    srand(time(nullptr)); // for a guaranteed connected graph
    for (int v = 2; v <= numNodes; v++) {
        int w = randomInt(1, v - 1); // Choose a random vertex to connect to
        while (w == v) { w = randomInt(1, numNodes); }
        graph[v][w] = 1;
        graph[w][v] = 1;
    }

    for (int k = 0; k < all.size(); k++) { // add more random edges
        for (int l = 0; l < all.size(); l++) {
            if (k + 1 != all[l]) {
                int rand = randomInt(0, 1);
                if (rand == 1) {
                    graph[k + 1][all[l]] = 1;
                    graph[all[l]][k + 1] = 1;
                }
            }
        }
    }
}

void generateSpecificGraph(int nodes, int edges) {
    srand(time(nullptr)); 

    // Generate random edges
    for (int i = 0; i < edges; ++i) {
        int u = randomInt(1, nodes); // Random source node
        int v = randomInt(1, nodes); // Random destination node

        // Avoid self-loops and duplicate edges
        while (u == v || graph[u][v] == 1) {
            u = randomInt(1, nodes);
            v = randomInt(1, nodes);
        }

        graph[u][v] = 1;
        graph[v][u] = 1;
    }
}

void printGraphEdges(vector<vector<int>> graph) {
    cout << "Graph Edges: { ";

    for (int i = 1; i <= N; i++) {
        for (int j = i + 1; j <= N; j++) {
            if (graph[i][j] == 1) {
                cout << "{" << i << ", " << j << "}, ";
            }
        }
    }

    cout << "}" << endl;
}

// Function to check if the given set of
// vertices in store array is a clique or not
bool is_clique(int b)
{
    vector<int> vertices;
    // Run a loop for all set of edges
    for (int i = 1; i < b; i++) {
        for (int j = i + 1; j < b; j++)

            // If any edge is missing
            if (graph[store[i]][store[j]] == 0) {
                vertices.clear();
                return false;
            }
            else {
                bool iExists = false;
                bool jExists = false;
                for (int f = 0; f < vertices.size(); f++) {
                    if (vertices[f] == store[i]) { iExists = true; }
                    if (vertices[f] == store[j]) { jExists = true; }
                }
                if (!iExists) { vertices.push_back(store[i]); }
                if (!jExists) { vertices.push_back(store[j]); }
            }
    }
    if (vertices.size() > BFclique.size()) {
        BFclique.clear();
        for (int k = 0; k < vertices.size(); k++) {
            BFclique.push_back(vertices[k]);
        }
    }
    return true;
}

int BFmaxCliques(int currIdx, int cliqueSize)
{
    // Maximal clique size
    int max_ = 0;

    // Check if any vertices from i+1 can be inserted
    for (int j = currIdx + 1; j <= N; j++) {

        // Add the vertex to store
        store[cliqueSize] = j;

        // If the graph is not a clique of size k then it cannot be a clique by adding another edge
        if (is_clique(cliqueSize + 1)) {

            // Update max
            max_ = max(max_, cliqueSize);

            // Check if another edge can be added
            max_ = max(max_, BFmaxCliques(j, cliqueSize + 1));
        }
    }
    return max_;
}

void HmaxCliques(vector<vector<int>> g) {

    //printGraphEdges(g);

    vector<int> vertex_degrees(g.size(), 0); //Find the degree of each vertex
    for (int i = 1; i < g.size(); i++) {
        for (int j = 1; j < g.size(); j++) {
            if (g[i][j] == 1) {
                vertex_degrees[i] += 1;
            }
        }
    }
    int num = 0;

    for (int f = 1; f < vertex_degrees.size(); f++) {
        if (vertex_degrees[f] != 0) { num++; }
    }

    bool isClique = true;
    for (int i = 1; i <= vertex_degrees.size() - 1; i++) { //Check if it forms a clique
        if (vertex_degrees[i] < num - 1 && vertex_degrees[i] != 0) {
            isClique = false;
            break;
        }
    }

    if (isClique) { // If it forms a clique, update the max clique value for each vertex if the new one is larger
        vector<int> v;

        for (int m = 0; m < vertex_degrees.size(); m++) { //store the vertices of the clique
            if (vertex_degrees[m] != 0) {
                v.push_back(m);
            }
        }

        if (v.size() > Hclique.size()) { // change the max clique vertices if the current one is larger
            Hclique.clear();
            for (int m = 0; m < v.size(); m++) {
                Hclique.push_back(v[m]);
            }
        }

        for (int i = 1; i <= vertex_degrees.size() - 1; i++) {
            if (vertex_degrees[i] != 0 && vertex_degrees[i] > maxC[i]) {
                maxC[i] = vertex_degrees[i] + 1;
            }
        }
    }

    else { // If not a clique
        //cout << "NOT a clique\n";
        int min = N;
        int minnode = 0;

        for (int i = 1; i < vertex_degrees.size(); i++) { //Find the minimum vertex
            if ((vertex_degrees[i] != 0) && (vertex_degrees[i] < min)) {
                min = vertex_degrees[i]; //ith vertex is the minimum
                minnode = i;
            }
        }


        vector<int> adjacent; //Find the index of vertices adjacent to min
        vector<vector<int>> subgraph = g;
        for (int i = 1; i < g.size(); i++) {
            if (g[minnode][i] == 1) { adjacent.push_back(i); }
        }


        for (int i = 1; i < N + 1; i++) { //create the subgraph with the adjacent nodes and min
            for (int j = 1; j < N + 1; j++) {
                bool isEdge = false;
                for (int f = 0; f < adjacent.size(); f++) {
                    if ((i == minnode && j == adjacent[f]) || (j == minnode && i == adjacent[f])) {
                        isEdge = true;
                    }
                }
                if (!isEdge) {
                    subgraph[i][j] = 0;
                }
            }
        }
        //cout << "called for subgraph\n";
        HmaxCliques(subgraph); // Check the subgraph

        for (int i = 1; i < g.size(); i++) { //Delete the minimum node from the graph (pruning)
            for (int j = 1; j < g.size(); j++) {
                if (i == minnode || j == minnode) {
                    g[i][j] = 0;
                }
            }
        }
        //cout << "called for pruned graph\n";
        if (!g.empty()) { // If the graph is not empty after pruning, keep going
            HmaxCliques(g);
        }
    }
}


int main() {
    srand(time(0));

    for (int i = 0; i < N + 1; i++) {

    }

    for (int i = 0; i <= N; i++) {
        for (int j = 0; j <= N; j++) {
            graph[i][j] = 0;
        }
    }

    generateRandomGraph(N);
    //generateSpecificGraph(N, 50);

    //initializing store

    for (int i = 0; i < N + 1; i++) {
        store[i] = 0;
    }

    //initializing d to keep track of vertex degrees

    for (int i = 1; i <= N; i++) {
        d[i] = 0;
    }

    for (int i = 1; i <= N; i++) {
        for (int j = 1; j <= N; j++) {
            if (graph[i][j] == 1) {
                d[i]++;
            }
        }
    }

    printGraphEdges(graph);

    auto start = std::chrono::high_resolution_clock::now();

    int maxSize = BFmaxCliques(0, 1);

    auto end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    cout << "Max Clique Size for Brute-Force Algorithm: " << maxSize << endl;
    cout << "Maximum Clique Vertices: \n";
    for (int i = 0; i < BFclique.size(); i++) {
        cout << "v" << BFclique[i] << " ";
    }
    cout << endl;

    cout << "Time taken for Brute Force Algorithm: " << duration.count() << " microseconds\n\n";

    auto start2 = std::chrono::high_resolution_clock::now();

    HmaxCliques(graph);

    auto end2 = std::chrono::high_resolution_clock::now();

    auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(end2 - start2);

    int hmax = 0;
    for (int i = 1; i <= N + 1; i++) {
        if (maxC[i] > hmax) {
            hmax = maxC[i];
        }
    }

    cout << "Max Clique Size for Heuristic Algorithm: " << hmax << endl;

    cout << "Maximum Clique Vertices: \n";
    for (int i = 0; i < Hclique.size(); i++) {
        cout << "v" << Hclique[i] << " ";
    }
    cout << endl;
    cout << "Time taken for Heuristic Algorithm: " << duration2.count() << " microseconds\n";

    return 0;
}
