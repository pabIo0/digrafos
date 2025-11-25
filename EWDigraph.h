#include <iostream>
#include <vector>
#include <list>
#include <stack>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm> 

using namespace std;

// Fila de Prioridade Minima Indexada
template <typename Key>
class IndexMinPQ {
private:
    int maxN;
    int n;
    vector<int> pq;     // heap binario usando indexaçao baseada em 1
    vector<int> qp;     // inverso de pq: qp[pq[i]] = pq[qp[i]] = i
    vector<Key> keys;   // chaves[i] = prioridade de i

    void exch(int i, int j) { swap(pq[i], pq[j]); qp[pq[i]] = i; qp[pq[j]] = j; }
    bool greater(int i, int j) const { return keys[pq[i]] > keys[pq[j]]; }
    void swim(int k) { while (k > 1 && greater(k/2, k)) { exch(k, k/2); k = k/2; } }
    void sink(int k) {
        while (2*k <= n) {
            int j = 2*k;
            if (j < n && greater(j, j+1)) j++;
            if (!greater(k, j)) break;
            exch(k, j);
            k = j;
        }
    }
    void validateIndex(int i) const {
        if (i < 0) throw out_of_range("indice e negativo");
        if (i >= maxN) throw out_of_range("indice >= capacidade maxima");
    }

public:
    IndexMinPQ(int N) : maxN(N), n(0), pq(N + 1), qp(N + 1, -1), keys(N + 1) {}

    bool isEmpty() const { return n == 0; }
    bool contains(int i) const { validateIndex(i); return qp[i] != -1; }

    void insert(int i, Key key) {
        validateIndex(i);
        if (contains(i)) throw logic_error("indice ja esta na fila de prioridade");
        n++;
        qp[i] = n;
        pq[n] = i;
        keys[i] = key;
        swim(n);
    }

    int delMin() {
        if (n == 0) throw out_of_range("Fila de prioridade vazia");
        int min = pq[1];
        exch(1, n--);
        sink(1);
        qp[min] = -1;
        return min;
    }

    void decreaseKey(int i, Key key) {
        validateIndex(i);
        if (!contains(i)) throw out_of_range("indice nao esta na fila de prioridade");
        if (keys[i] < key)
            throw invalid_argument("decreaseKey() com chave maior que a atual");
        keys[i] = key;
        swim(qp[i]);
    }
};

/**
 * @struct Edge
 * @brief Representa uma aresta direcionada e ponderada (v -> w).
 */
struct Edge {
    int v, w;
    double weight;

    Edge(int v = -1, int w = -1, double weight = 0.0) : v(v), w(w), weight(weight) {}

    int from() const { return v; }
    int to() const { return w; }
};

/**
 * @class EWDigraph
 * @brief Representa um digrafo ponderado usando lista de adjacência.
 */
class EWDigraph {
private:
    int V;
    int E;
    vector<list<Edge>> adj;

public:
    // Construtor a partir do numero de vertices
    EWDigraph(int V) : V(V), E(0) {
        if (V < 0) throw invalid_argument("Numero de vertices nao pode ser negativo");
        adj.resize(V);
    }

    // Construtor a partir de um arquivo de entrada
    EWDigraph(istream &in) {
        if (!in) throw invalid_argument("Stream de entrada nulo.");
        in >> V;
        adj.resize(V);
        E = 0;

        int totalEdges;
        in >> totalEdges;

        for (int i = 0; i < totalEdges; i++) {
            int v, w;
            double weight;
            in >> v >> w >> weight;
            addEdge(Edge(v, w, weight));
        }
    }

    // Adiciona uma aresta direcionada v -> w
    void addEdge(Edge e) {
        E++;
        adj[e.from()].push_back(e);
    }

    int getV() const { return V; }
    int getE() const { return E; }


     // Retorna todas as arestas do grafo 
    vector<Edge> getAllEdges() const {
        vector<Edge> edges;
        for (int v = 0; v < V; ++v) {
            for (const auto& edge : adj[v]) {
                edges.push_back(edge);
            }
        }
        return edges;
    }
     

     list<Edge> getAdj(int v) const
    {
        return adj[v];
    }

    //Mostra o Grafo
    void show(){
       vector<Edge> edges;
       cout << "Grafo com " << V << " vertices e " << E << " arestas." << endl;
       cout << "---------------------------------------------" << endl;
       for (int v = 0; v < V; ++v) {
           for (const auto& edge : adj[v]) {
                cout << "  " << v << " -> " << edge.w << ": "<< edge.weight << "\n";
           }
       }
        cout << "---------------------------------------------" << endl;
    }

    //Mostra o grafo no formato dot
    void showDot(){
        cout << "digraph {\n";
        cout << "  node [shape=circle];\n";
        cout << "  edge [labeldistance=1.5];\n";
        for (int v = 0; v < V; ++v) {
            for (const auto& edge : adj[v]) {
                cout << "  " << v << " -> " << edge.w
                        << " [label=\"" << fixed << setprecision(2) << edge.weight << "\"];\n";
            }
        }
            cout << "}\n";
    }

    // Iterador para percorrer as arestas adjacentes a um vertice
    class adjIterator {
    private:
        const EWDigraph &G;
        typename list<Edge>::const_iterator it;
    public:
        adjIterator(const EWDigraph &G, int v) : G(G) {
            it = G.adj[v].begin();
        }

        Edge beg(int v) {
            it = G.adj[v].begin();
            return (it != G.adj[v].end()) ? *it : Edge(-1, -1);
        }

        Edge nxt(int v) {
            if (it != G.adj[v].end()) ++it;
            return (it != G.adj[v].end()) ? *it : Edge(-1, -1);
        }

        bool end(int v) const {
            return it == G.adj[v].end();
        }
    };
};


