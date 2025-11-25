#include <iostream>
#include <vector>
#include <list>
#include <stack>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <limits>
#include <cmath>

using namespace std;

// --- Fila de Prioridade Mínima Indexada (Necessária para Dijkstra) ---
template <typename Key>
class IndexMinPQ {
private:
    int maxN;
    int n;
    vector<int> pq;     // heap binário base 1
    vector<int> qp;     // inverso de pq
    vector<Key> keys;   // prioridades

    void exch(int i, int j) { swap(pq[i], pq[j]); qp[pq[i]] = i; qp[pq[j]] = j; }
    bool greater(int i, int j) const { return keys[pq[i]] > keys[pq[j]]; }
    
    void swim(int k) { 
        while (k > 1 && greater(k/2, k)) { exch(k, k/2); k = k/2; } 
    }
    
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
        if (i < 0) throw out_of_range("Indice negativo");
        if (i >= maxN) throw out_of_range("Indice maior que a capacidade");
    }

public:
    IndexMinPQ(int N) : maxN(N), n(0), pq(N + 1), qp(N + 1, -1), keys(N + 1) {}

    bool isEmpty() const { return n == 0; }
    bool contains(int i) const { validateIndex(i); return qp[i] != -1; }

    void insert(int i, Key key) {
        validateIndex(i);
        if (contains(i)) throw logic_error("Indice ja existe na fila");
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
        if (!contains(i)) throw out_of_range("Indice nao existe na fila");
        if (keys[i] < key) throw invalid_argument("Chave nova eh maior que a atual");
        keys[i] = key;
        swim(qp[i]);
    }
};

// --- Estrutura da Aresta (Atualizada para Fluxo) ---
struct Edge {
    int v, w;
    double weight; // Capacidade/Custo
    double flow;   // Fluxo atual (Novo!)

    Edge(int v = -1, int w = -1, double weight = 0.0) 
        : v(v), w(w), weight(weight), flow(0.0) {}

    int from() const { return v; }
    int to() const { return w; }

    // Retorna o vértice oposto (essencial para grafos residuais)
    int other(int vertex) const {
        if (vertex == v) return w;
        if (vertex == w) return v;
        throw invalid_argument("Vertice invalido na aresta");
    }

    // Métodos auxiliares para Ford-Fulkerson
    double capacity() const { return weight; }
    
    double residualCapacityTo(int vertex) const {
        if (vertex == w) return weight - flow; // Forward: capacidade restante
        else if (vertex == v) return flow;     // Backward: fluxo que pode retornar
        else throw invalid_argument("Vertice invalido na aresta");
    }

    void addResidualFlowTo(int vertex, double delta) {
        if (vertex == w) flow += delta;       // Aumenta fluxo
        else if (vertex == v) flow -= delta;  // Diminui fluxo
        else throw invalid_argument("Vertice invalido na aresta");
    }
};

// --- Classe do Dígrafo Ponderado ---
class EWDigraph {
private:
    int V;
    int E;
    vector<list<Edge>> adj;

public:
    EWDigraph(int V) : V(V), E(0) {
        if (V < 0) throw invalid_argument("Numero de vertices nao pode ser negativo");
        adj.resize(V);
    }

    EWDigraph(istream &in) {
        if (!in) throw invalid_argument("Stream de entrada invalido");
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

    void addEdge(Edge e) {
        E++;
        adj[e.from()].push_back(e);
    }

    int getV() const { return V; }
    int getE() const { return E; }

    // Retorna todas as arestas (útil para clonar)
    vector<Edge> getAllEdges() const {
        vector<Edge> edges;
        for (int v = 0; v < V; ++v) {
            for (const auto& edge : adj[v]) {
                edges.push_back(edge);
            }
        }
        return edges;
    }

    // Retorna referência para permitir modificação do fluxo
    list<Edge>& getAdj(int v) { return adj[v]; }
    const list<Edge>& getAdj(int v) const { return adj[v]; }

    // Verifica se há arestas negativas (para proteção do Dijkstra)
    bool hasNegativeEdge() const {
        for (int v = 0; v < V; ++v) {
            for (const auto& e : adj[v]) {
                if (e.weight < 0) return true;
            }
        }
        return false;
    }

    void show() {
        cout << "Grafo com " << V << " vertices e " << E << " arestas." << endl;
        cout << "---------------------------------------------" << endl;
        for (int v = 0; v < V; ++v) {
            for (const auto& edge : adj[v]) {
                if (abs(edge.weight) > 1e-9 || edge.flow > 0)
                    cout << "  " << v << " -> " << edge.w << ": " << edge.weight << "\n";
            }
        }
        cout << "---------------------------------------------" << endl;
    }
    
    void showDot() {
        cout << "digraph {\n";
        cout << "  node [shape=circle];\n";
        cout << "  edge [labeldistance=1.5];\n";
        for (int v = 0; v < V; ++v) {
            for (const auto& edge : adj[v]) {
                if (abs(edge.weight) > 1e-9) {
                    cout << "  " << v << " -> " << edge.w
                         << " [label=\"" << fixed << setprecision(2) << edge.weight << "\"];\n";
                }
            }
        }
        cout << "}\n";
    }
};