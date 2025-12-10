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
#include <cmath> // Necessário para abs()

using namespace std;

// =================================================================
// ESTRUTURA DE DADOS AUXILIAR: Fila de Prioridade Mínima Indexada
// =================================================================
// OBJETIVO: Permitir que o algoritmo de DIJKSTRA rode em tempo O(E log V).
// DIFERENÇA PARA std::priority_queue: Esta estrutura permite alterar a 
// prioridade de um item que já está na fila (decreaseKey), algo que a
// fila padrão do C++ não suporta nativamente de forma eficiente.
//
// ESTRUTURA: Heap Binário implementado sobre vetores.
// - pq[]: O Heap em si (base 1). Guarda os índices dos vértices.
// - qp[]: Mapeamento inverso. Diz onde o vértice 'i' está dentro do heap 'pq'.
// - keys[]: Guarda os valores de prioridade (as distâncias distTo).
// =================================================================
template <typename Key>
class IndexMinPQ {
private:
    int maxN;           // Capacidade máxima de elementos
    int n;              // Número atual de elementos na fila
    vector<int> pq;     // Heap Binário (indices)
    vector<int> qp;     // Array Inverso: qp[i] = posição de i no heap pq
    vector<Key> keys;   // Prioridades (chaves)

    // Auxiliar: Troca dois elementos no heap e atualiza o mapeamento inverso
    void exch(int i, int j) { swap(pq[i], pq[j]); qp[pq[i]] = i; qp[pq[j]] = j; }
    
    // Auxiliar: Compara a prioridade de dois índices no heap
    bool greater(int i, int j) const { return keys[pq[i]] > keys[pq[j]]; }
    
    // Auxiliar: SWIM (Nadar) - Corrige a ordem do heap de baixo para cima.
    // Usado quando inserimos um novo elemento ou diminuímos uma chave.
    void swim(int k) { 
        while (k > 1 && greater(k/2, k)) { exch(k, k/2); k = k/2; } 
    }
    
    // Auxiliar: SINK (Afundar) - Corrige a ordem do heap de cima para baixo.
    // Usado quando removemos o topo (menor elemento).
    void sink(int k) {
        while (2*k <= n) {
            int j = 2*k;
            if (j < n && greater(j, j+1)) j++; // Escolhe o filho menor
            if (!greater(k, j)) break;         // Se a ordem está certa, para
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
    
    // Verifica se o índice i está na fila em tempo O(1) usando o array qp
    bool contains(int i) const { validateIndex(i); return qp[i] != -1; }

    // Insere um índice com uma prioridade associada
    void insert(int i, Key key) {
        validateIndex(i);
        if (contains(i)) throw logic_error("Indice ja existe na fila");
        n++;
        qp[i] = n;
        pq[n] = i;
        keys[i] = key;
        swim(n); // Restaura a propriedade do heap subindo o elemento
    }

    // Remove e retorna o índice de menor prioridade (topo do heap)
    int delMin() {
        if (n == 0) throw out_of_range("Fila de prioridade vazia");
        int min = pq[1];
        exch(1, n--); // Troca o topo com o último e reduz tamanho
        sink(1);      // Restaura a propriedade do heap descendo o novo topo
        qp[min] = -1; // Marca como removido
        return min;
    }

    // Atualiza o valor da chave associada ao índice i
    // CRUCIAL PARA DIJKSTRA: Quando encontramos um caminho mais curto para 'i',
    // chamamos isso para atualizar sua posição na fila.
    void decreaseKey(int i, Key key) {
        validateIndex(i);
        if (!contains(i)) throw out_of_range("Indice nao existe na fila");
        if (keys[i] < key) throw invalid_argument("Chave nova eh maior que a atual");
        keys[i] = key;
        swim(qp[i]); // A prioridade diminuiu, então o elemento deve "subir" no heap
    }
};

// =================================================================
// ESTRUTURA DA ARESTA (Híbrida: Caminho Mínimo + Fluxo)
// =================================================================
struct Edge {
    int v, w;      // Aresta direcionada de v -> w
    double weight; // CAMINHOS MÍNIMOS: Representa o Custo/Peso.
                   // FLUXO MÁXIMO: Representa a CAPACIDADE da aresta.
    
    double flow;   // FLUXO MÁXIMO: Representa o fluxo atual passando pela aresta.
                   // (Não é usado nos algoritmos de caminho mínimo).

    Edge(int v = -1, int w = -1, double weight = 0.0) 
        : v(v), w(w), weight(weight), flow(0.0) {}

    int from() const { return v; }
    int to() const { return w; }

    // --- MÉTODOS PARA NAVEGAÇÃO EM GRAFOS (ESSENCIAL PARA GRAFO RESIDUAL) ---
    
    // Retorna o vértice na outra ponta da aresta, dado um vértice conhecido.
    // Útil quando navegamos "contra a mão" no grafo residual.
    int other(int vertex) const {
        if (vertex == v) return w;
        if (vertex == w) return v;
        throw invalid_argument("Vertice invalido na aresta");
    }

    // --- MÉTODOS AUXILIARES PARA FORD-FULKERSON ---

    // Retorna a capacidade total (semântica de fluxo)
    double capacity() const { return weight; }
    
    // Calcula a Capacidade Residual em direção a um vértice 'vertex'.
    // Lógica do Grafo Residual:
    // 1. Se vamos na direção original (v->w): Capacidade restante = (Capacidade - Fluxo Atual)
    // 2. Se vamos na direção oposta (w->v): Capacidade de retorno = (Fluxo Atual)
    //    (Isso representa a capacidade de "desfazer" ou redirecionar fluxo já enviado)
    double residualCapacityTo(int vertex) const {
        if (vertex == w) return weight - flow; // Forward edge
        else if (vertex == v) return flow;     // Backward edge
        else throw invalid_argument("Vertice invalido na aresta");
    }

    // Atualiza o fluxo ao passar por esta aresta no grafo residual
    // delta: quantidade de fluxo sendo empurrada pelo caminho aumentante
    void addResidualFlowTo(int vertex, double delta) {
        if (vertex == w) flow += delta;       // Forward: aumenta o fluxo usado
        else if (vertex == v) flow -= delta;  // Backward: diminui o fluxo usado (cancela fluxo)
        else throw invalid_argument("Vertice invalido na aresta");
    }
};

// =================================================================
// CLASSE DO DÍGRAFO PONDERADO
// =================================================================
class EWDigraph {
private:
    int V; // Número de vértices
    int E; // Número de arestas
    vector<list<Edge>> adj; // Lista de Adjacência: adj[v] guarda as arestas saindo de v

public:
    // Construtor básico
    EWDigraph(int V) : V(V), E(0) {
        if (V < 0) throw invalid_argument("Numero de vertices nao pode ser negativo");
        adj.resize(V);
    }

    // Construtor que lê de arquivo
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

    // Retorna todas as arestas do grafo (útil para clonagem ou iteração global)
    vector<Edge> getAllEdges() const {
        vector<Edge> edges;
        for (int v = 0; v < V; ++v) {
            for (const auto& edge : adj[v]) {
                edges.push_back(edge);
            }
        }
        return edges;
    }

    // Acesso à lista de adjacência (referência para permitir alterar fluxo das arestas)
    list<Edge>& getAdj(int v) { return adj[v]; }
    const list<Edge>& getAdj(int v) const { return adj[v]; }

    // Auxiliar: Verifica segurança para rodar Dijkstra
    // Dijkstra falha (loop ou resultado errado) com arestas negativas.
    bool hasNegativeEdge() const {
        for (int v = 0; v < V; ++v) {
            for (const auto& e : adj[v]) {
                if (e.weight < 0) return true;
            }
        }
        return false;
    }

    // Exibe o grafo no terminal
    void show() {
        cout << "Grafo com " << V << " vertices e " << E << " arestas." << endl;
        cout << "---------------------------------------------" << endl;
        for (int v = 0; v < V; ++v) {
            for (const auto& edge : adj[v]) {
                // Filtro Visual: 
                // Mostra a aresta se ela tiver peso significativo (original) OU fluxo (usada no algoritmo)
                // 'abs(edge.weight) > 1e-9': Garante que mostramos arestas negativas, mas escondemos
                // as arestas virtuais reversas de capacidade 0 criadas pelo Ford-Fulkerson.
                if (abs(edge.weight) > 1e-9 || edge.flow > 0)
                    cout << "  " << v << " -> " << edge.w << ": " << edge.weight << "\n";
            }
        }
        cout << "---------------------------------------------" << endl;
    }
    
    // Graphviz
    void showDot() {
        cout << "digraph {\n";
        cout << "  node [shape=circle];\n";
        cout << "  edge [labeldistance=1.5];\n";
        for (int v = 0; v < V; ++v) {
            for (const auto& edge : adj[v]) {
                // Filtra arestas artificiais de capacidade 0
                if (abs(edge.weight) > 1e-9) {
                    cout << "  " << v << " -> " << edge.w
                         << " [label=\"" << fixed << setprecision(2) << edge.weight << "\"];\n";
                }
            }
        }
        cout << "}\n";
    }
};