#include "EWDigraph.h"
#include <functional>
#include <queue>
#include <stack>
#include <iomanip>
#include <limits>

using namespace std;

class Process
{
private:
    EWDigraph G;

public:
    Process(const EWDigraph &G) : G(G) {}

    // ==========================================================
    // 1. ALGORITMO DE DIJKSTRA
    // ----------------------------------------------------------
    // TIPO: Algoritmo Guloso (Greedy).
    // USO: Caminhos mínimos com pesos NÃO-NEGATIVOS.
    // COMPLEXIDADE: O(E log V) com Fila de Prioridade Indexada.
    // ==========================================================
    class DjkstraSP
    {
    private:
        vector<double> distTo; // Menor distância conhecida da origem 's' até 'v'
        vector<Edge> edgeTo;   // Última aresta usada para chegar em 'v' (para reconstruir o caminho)
        IndexMinPQ<double> pq; // Fila de prioridade para escolher o próximo vértice mais próximo
      
    public:
        DjkstraSP(const EWDigraph &G, int s) 
            : distTo(G.getV(), numeric_limits<double>::infinity()), 
              edgeTo(G.getV()), 
              pq(G.getV())
        {
            distTo[s] = 0.0;    
            pq.insert(s, 0.0);
            
            // Enquanto houver vértices na fila (fronteira de exploração)
            while (!pq.isEmpty())
            {
                // Remove o vértice com a menor distância estimada (Guloso)
                int v = pq.delMin();
                
                // Relaxa todas as arestas que saem deste vértice
                for (const Edge &e : G.getAdj(v))
                {
                    relax(e);
                }
            }
        }   

        // RELAXAMENTO DE ARESTA (Conceito Chave)
        // Pergunta: "Se eu for de s->v->w, é mais rápido do que o caminho atual s->w?"
        void relax(const Edge &e)
        {
            int v = e.from();
            int w = e.to();
            // Se o caminho passando por 'v' for menor que o que eu já conheço para 'w'...
            if (distTo[w] > distTo[v] + e.weight)
            {
                distTo[w] = distTo[v] + e.weight; // Atualiza a distância
                edgeTo[w] = e;                    // Atualiza o "pai" na árvore de caminhos
                
                // Atualiza a prioridade na fila ou insere se não existir
                if (pq.contains(w)) pq.decreaseKey(w, distTo[w]);
                else pq.insert(w, distTo[w]);
            }
        }

        double distanceTo(int v) const { return distTo[v]; }
        bool hasPathTo(int v) const { return distTo[v] < numeric_limits<double>::infinity(); }
        
        // Reconstrói o caminho empilhando as arestas de edgeTo[] de trás para frente
        stack<Edge> pathTo(int v) const
        {   
            stack<Edge> path;
            if (!hasPathTo(v)) return path;
            for (Edge e = edgeTo[v]; e.from() != -1; e = edgeTo[e.from()])
                path.push(e);
            return path;
        }       
    };

    // ==========================================================
    // 2. ACYCLIC SP (Caminhos Mínimos em DAGs)
    // ----------------------------------------------------------
    // TIPO: Programação Dinâmica / Ordenação Topológica.
    // USO: Apenas Grafos ACÍCLICOS (DAG). Aceita pesos negativos.
    // COMPLEXIDADE: O(V + E) - Linear (muito rápido).
    // ==========================================================
    class AcyclicSP
    {
    private:
        vector<double> distTo;
        vector<Edge> edgeTo;
        bool cycleFound; // Flag de segurança (Fail-Fast)

    public:
        AcyclicSP(const EWDigraph &G, int s) : cycleFound(false)
        {
            distTo.resize(G.getV(), numeric_limits<double>::infinity());
            edgeTo.resize(G.getV(), Edge(-1, -1, 0.0));
            distTo[s] = 0.0;

            // Variáveis para a DFS (Busca em Profundidade)
            vector<bool> marked(G.getV(), false);
            vector<bool> onStack(G.getV(), false); // Para detectar back-edges (ciclos)
            stack<int> reversePost; // Armazena a Ordem Topológica Inversa
            
            // Função Lambda para DFS com Detecção de Ciclo
            function<void(int)> dfs = [&](int v) {
                marked[v] = true;
                onStack[v] = true; // Entrando na pilha de recursão
                
                for (const Edge &e : G.getAdj(v)) {
                    if (cycleFound) return;
                    int w = e.to();
                    if (!marked[w]) {
                        dfs(w);
                    } else if (onStack[w]) {
                        // Se encontramos um vizinho que já está na pilha de recursão atual,
                        // significa que há um caminho de volta -> CICLO!
                        cycleFound = true; 
                    }
                }
                onStack[v] = false; // Saindo da pilha de recursão
                reversePost.push(v); // Adiciona na pilha (Pós-ordem reversa = Ordem Topológica)
            };

            // Executa DFS em todos os componentes
            for (int v = 0; v < G.getV(); v++)
                if (!marked[v] && !cycleFound) dfs(v);

            // SE HOUVER CICLO, O ALGORITMO PARA AQUI.
            // AcyclicSP não funciona com ciclos.
            if (cycleFound) return;

            // Processa vértices na Ordem Topológica
            // Garante que quando processamos 'v', todas as arestas chegando em 'v' já foram relaxadas.
            while (!reversePost.empty())
            {
                int v = reversePost.top();
                reversePost.pop();
                for (const Edge &e : G.getAdj(v))
                    relax(e);
            }
        }

        void relax(const Edge &e)
        {
            int v = e.from();
            int w = e.to();
            if (distTo[w] > distTo[v] + e.weight)
            {
                distTo[w] = distTo[v] + e.weight;
                edgeTo[w] = e;
            }
        }

        bool isValid() const { return !cycleFound; }
        double distanceTo(int v) const { return distTo[v]; }
        bool hasPathTo(int v) const { return distTo[v] < numeric_limits<double>::infinity(); }
        
        stack<Edge> pathTo(int v) const
        {   
            stack<Edge> path;
            if (!hasPathTo(v)) return path;
            for (Edge e = edgeTo[v]; e.from() != -1; e = edgeTo[e.from()])
                path.push(e);
            return path;
        }
    };

    // ==========================================================
    // 3. BELLMAN-FORD (Algoritmo Robusto)
    // ----------------------------------------------------------
    // TIPO: Programação Dinâmica.
    // USO: Pesos negativos permitidos. Detecta CICLOS NEGATIVOS.
    // COMPLEXIDADE: O(V * E) - Mais lento que Dijkstra.
    // ==========================================================
    class BellmanFordSP
    {
    private:
        vector<double> distTo;
        vector<Edge> edgeTo;
        stack<Edge> cycle; // Armazena o ciclo negativo se encontrado
        bool hasCycle;

    public:
        BellmanFordSP(const EWDigraph &G, int s) : hasCycle(false)
        {
            distTo.resize(G.getV(), numeric_limits<double>::infinity());
            edgeTo.resize(G.getV(), Edge(-1, -1, 0.0));
            distTo[s] = 0.0;    

            // PASSO 1: Relaxar todas as arestas |V|-1 vezes.
            // Teoria: O caminho mínimo simples mais longo tem no máximo |V|-1 arestas.
            for (int i = 0; i < G.getV() - 1; i++)
            {
                for (int v = 0; v < G.getV(); v++)
                {
                    for (const Edge &e : G.getAdj(v))
                    {
                        relax(e);
                    }
                }
            }

            // PASSO 2: Verificar Ciclos Negativos (V-ésima iteração).
            // Se após V-1 iterações ainda for possível diminuir alguma distância,
            // então matematicamente existe um ciclo negativo acessível.
            for (int v = 0; v < G.getV(); v++)
            {
                for (const Edge &e : G.getAdj(v))
                {
                    int w = e.to();
                    // Usamos uma tolerância (+ 1e-9) para evitar erros de ponto flutuante (double)
                    // onde 0.0000000001 seria considerado uma "melhora".
                    if (distTo[v] != numeric_limits<double>::infinity() && 
                        distTo[w] > distTo[v] + e.weight + 1e-9)
                    {
                        hasCycle = true;
                        
                        // LÓGICA PARA RECONSTRUIR O CICLO NEGATIVO:
                        // 1. Voltamos V vezes no edgeTo[] para garantir que entramos no loop do ciclo.
                        int cycleNode = w;
                        for(int i = 0; i < G.getV(); ++i) {
                           if (edgeTo[cycleNode].from() == -1) break; 
                           cycleNode = edgeTo[cycleNode].from();
                        }

                        // 2. Percorremos o edgeTo[] até encontrar o mesmo vértice novamente.
                        stack<Edge> tempCycle;
                        int v_curr = cycleNode;
                        while (true) {
                            Edge edge = edgeTo[v_curr];
                            if (edge.from() == -1) break;
                            tempCycle.push(edge);
                            v_curr = edge.from();
                            if (v_curr == cycleNode && !tempCycle.empty()) {
                                cycle = tempCycle;
                                return;
                            }
                        }
                        return;
                    }
                }
            }
        }   

        void relax(const Edge &e)
        {
            int v = e.from();
            int w = e.to();
            if (distTo[v] != numeric_limits<double>::infinity() && 
                distTo[w] > distTo[v] + e.weight)
            {
                distTo[w] = distTo[v] + e.weight;
                edgeTo[w] = e;
            }
        }

        bool hasNegativeCycle() const { return hasCycle; }
        stack<Edge> getNegativeCycle() const { return cycle; } // Retorna a pilha com o ciclo
        
        double distanceTo(int v) const { return distTo[v]; }
        bool hasPathTo(int v) const { return distTo[v] < numeric_limits<double>::infinity(); }
        
        stack<Edge> pathTo(int v) const
        {   
            stack<Edge> path;
            if (!hasPathTo(v)) return path;
            for (Edge e = edgeTo[v]; e.from() != -1; e = edgeTo[e.from()])
                path.push(e);
            return path;
        }
    };

    // ==========================================================
    // 4. FLOYD-WARSHALL (All-Pairs Shortest Path)
    // ----------------------------------------------------------
    // TIPO: Programação Dinâmica.
    // USO: Menor caminho entre TODOS os pares de vértices.
    // COMPLEXIDADE: O(V^3) - Cúbica.
    // ==========================================================
    class FloydWarshall
    {
    private:
        vector<vector<double>> dist; // Matriz de Distâncias
        bool hasNegativeCycle;

    public:
        FloydWarshall(const EWDigraph &G) : hasNegativeCycle(false)
        {
            int V = G.getV();
            dist.resize(V, vector<double>(V, numeric_limits<double>::infinity()));

            // Inicialização da Matriz
            // dist[i][i] = 0
            // dist[i][j] = peso da aresta se existir
            for (int v = 0; v < V; v++) {
                dist[v][v] = 0.0;
                for (const Edge &e : G.getAdj(v)) {
                    dist[e.from()][e.to()] = e.weight;
                }
            }

            // O Coração do Algoritmo: 3 loops aninhados
            // k: vértice intermediário permitido
            // i: origem
            // j: destino
            // Relação: dist[i][j] = min(dist[i][j], dist[i][k] + dist[k][j])
            for (int k = 0; k < V; k++) {
                for (int i = 0; i < V; i++) {
                    for (int j = 0; j < V; j++) {
                        if (dist[i][k] != numeric_limits<double>::infinity() &&
                            dist[k][j] != numeric_limits<double>::infinity()) {
                            
                            // Tolerância 1e-9 aplicada aqui também
                            if (dist[i][j] > dist[i][k] + dist[k][j] + 1e-9) {
                                dist[i][j] = dist[i][k] + dist[k][j];
                            }
                        }
                    }
                }
                // Se a distância de k para k for negativa, há um ciclo negativo passando por k
                if (dist[k][k] < -1e-9) { hasNegativeCycle = true; return; }
            }
        }

        double distTo(int s, int t) { return dist[s][t]; }
        bool hasNegCycle() { return hasNegativeCycle; }
    };

    // ==========================================================
    // 5. FORD-FULKERSON (Fluxo Máximo / Min-Cut)
    // ----------------------------------------------------------
    // TIPO: Método de Caminhos Aumentantes (Edmonds-Karp usa BFS).
    // USO: Calcular fluxo máximo em uma rede.
    // LÓGICA: Encontrar caminhos no "Grafo Residual" e saturar arestas.
    // ==========================================================
    class FordFulkerson
    {
    private:
        vector<bool> marked;  // Visitados na BFS (ajuda a identificar o Min-Cut)
        vector<Edge*> edgeTo; // Caminho encontrado pela BFS
        double value;         // Valor total do fluxo

        // BFS no Grafo Residual
        // Retorna true se existe um caminho de s -> t onde capacidade residual > 0
        bool hasAugmentingPath(EWDigraph &G, int s, int t) {
            edgeTo.assign(G.getV(), nullptr);
            marked.assign(G.getV(), false);
            queue<int> q;
            q.push(s);
            marked[s] = true;

            while (!q.empty() && !marked[t]) {
                int v = q.front(); q.pop();
                for (Edge &e : G.getAdj(v)) {
                    int w = e.other(v); // Pega o vértice oposto (seja forward ou backward)
                    
                    // Se há capacidade residual na aresta E o vizinho não foi visitado
                    if (e.residualCapacityTo(w) > 0 && !marked[w]) {
                        edgeTo[w] = &e; // Guarda a aresta usada
                        marked[w] = true;
                        q.push(w);
                    }
                }
            }
            return marked[t]; // Chegamos no destino?
        }

    public:
        FordFulkerson(EWDigraph &G, int s, int t) : value(0.0) {
            // PRÉ-PROCESSAMENTO:
            // O algoritmo precisa navegar "para trás" nas arestas.
            // O EWDigraph original só tem arestas de ida.
            // Aqui, adicionamos as arestas reversas virtuais (com capacidade 0) se não existirem.
            vector<Edge> originalEdges = G.getAllEdges();
            for(const auto& e : originalEdges) {
                bool hasReverse = false;
                for(const auto& rev : G.getAdj(e.to())) {
                    if(rev.to() == e.from()) { hasReverse = true; break; }
                }
                if(!hasReverse) {
                    // Adiciona aresta reversa (peso/capacidade = 0)
                    G.addEdge(Edge(e.to(), e.from(), 0.0));
                }
            }

            // ENQUANTO EXISTIR CAMINHO AUMENTANTE (BFS):
            while (hasAugmentingPath(G, s, t)) {
                // 1. Calcular o Gargalo (Bottleneck) do caminho encontrado
                double bottle = numeric_limits<double>::infinity();
                for (int v = t; v != s; v = edgeTo[v]->other(v)) {
                    bottle = min(bottle, edgeTo[v]->residualCapacityTo(v));
                }

                // 2. Atualizar o Fluxo Residual ao longo do caminho
                // (Adiciona fluxo na direção forward, remove na backward)
                for (int v = t; v != s; v = edgeTo[v]->other(v)) {
                    edgeTo[v]->addResidualFlowTo(v, bottle);
                }

                value += bottle; // Soma ao fluxo total
            }
        }

        double getValue() const { return value; }
        
        // Retorna true se o vértice 'v' está no conjunto S do Corte Mínimo (alcançável de s no grafo residual)
        bool inCut(int v) const { return marked[v]; }
    };
};

int main(int argc, char *argv[])
{
    if (argc < 2) {
        cerr << "Uso: " << argv[0] << " <arquivo_do_grafo>" << endl;
        return 1;
    }
    ifstream in(argv[1]);
    if (!in) {
        cerr << "Erro ao abrir o arquivo." << endl;
        return 1;
    }

    EWDigraph G(in);
    G.show();
    // G.showDot(); // Graphviz
    
    int s = 0; // Vértice de origem padrão

    // ==========================================================
    // 1. ALGORITMO DE DIJKSTRA
    // ----------------------------------------------------------
    // TIPO: Algoritmo Guloso (Greedy Algorithm).
    // USO: Encontrar caminhos mínimos de uma origem 's' para todos os outros vértices.
    //      **RESTRIÇÃO:** Só funciona se todas as arestas tiverem PESOS NÃO-NEGATIVOS.
    // COMPLEXIDADE: O(E log V) (usando Binary Heap / IndexMinPQ).
    //
    // PECULIARIDADES:
    // 1. Abordagem Gulosa: Assume que, ao fechar um vértice, já encontrou o melhor caminho
    //    possível para ele. Por isso falha com pesos negativos (não consegue "voltar atrás").
    // 2. Fila de Prioridade: Essencial para garantir a complexidade logarítmica, escolhendo
    //    sempre o vértice não visitado mais próximo da origem.
    // 3. Relaxamento: A operação chave é verificar se passar por 'v' melhora o caminho até 'w'.
    //    distTo[w] > distTo[v] + weight.
    // ----------------------------------------------------------
    cout << "\n------- Dijkstra SP (s=" << s << ") -------" << endl;
    // Proteção: Dijkstra falha com arestas negativas
    if (G.hasNegativeEdge()) {
        cout << "Aviso: Grafo contem arestas negativas. Dijkstra nao e confiavel." << endl;
    } else {
        Process::DjkstraSP dsp(G, s);
        for (int v = 0; v < G.getV(); v++) {
            if (dsp.hasPathTo(v)) {
                cout << "Caminho de " << s << " para " << v << " (" << fixed << setprecision(2) << dsp.distanceTo(v) << "): ";
                stack<Edge> path = dsp.pathTo(v);
                while(!path.empty()) {
                    Edge e = path.top(); path.pop();
                    cout << e.from() << "->" << e.to() << " (" << e.weight << ") ";
                }
                cout << endl;
            } else {
                cout << "Nao ha caminho de " << s << " para " << v << endl;
            }
        }
    }

    // ==========================================================
    // 2. ACYCLIC SP (Caminhos Mínimos em DAGs)
    // ----------------------------------------------------------
    // TIPO: Programação Dinâmica / Ordenação Topológica.
    // USO: Caminhos mínimos em Grafos Direcionados Acíclicos (DAGs).
    //      Funciona com PESOS NEGATIVOS, desde que não haja ciclos.
    // COMPLEXIDADE: O(V + E) - Tempo Linear (o mais rápido possível).
    //
    // PECULIARIDADES:
    // 1. Ordenação Topológica: O segredo é relaxar os vértices numa ordem específica
    //    (onde todas as dependências vêm antes). Isso garante que só precisamos visitar
    //    cada aresta uma única vez.
    // 2. Verificação de Ciclo: Obrigatória. Se o grafo tiver um ciclo, não existe ordem
    //    topológica e o algoritmo falha/produz lixo.
    // 3. DFS (Depth First Search): Usada internamente para gerar a ordem topológica.
    // ----------------------------------------------------------
    cout << "\n------- Acyclic SP (s=" << s << ") -------" << endl;
    Process::AcyclicSP acyclicSP(G, s);
    // Proteção: AcyclicSP falha com Ciclos
    if (acyclicSP.isValid()) {
        for (int v = 0; v < G.getV(); v++) {
            if (acyclicSP.hasPathTo(v))
                 cout << "Caminho de " << s << " para " << v << " (" << fixed << setprecision(2) << acyclicSP.distanceTo(v) << ")" << endl;
            else
                cout << "Nao ha caminho de " << s << " para " << v << endl;
        }
    } else {
        cout << "Erro: O grafo contem ciclos! AcyclicSP nao aplicavel." << endl;
    }

    // ==========================================================
    // 3. BELLMAN-FORD (Algoritmo Robusto)
    // ----------------------------------------------------------
    // TIPO: Programação Dinâmica (Relaxamento Iterativo).
    // USO: Caminhos mínimos em grafos gerais.
    //      **VANTAGEM:** Suporta PESOS NEGATIVOS e detecta CICLOS NEGATIVOS.
    // COMPLEXIDADE: O(V * E) - Muito mais lento que Dijkstra.
    //
    // PECULIARIDADES:
    // 1. Força Bruta Inteligente: Relaxa TODAS as arestas V vezes.
    // 2. Teorema do Caminho Simples: Um caminho mínimo sem ciclos tem no máximo V-1 arestas.
    //    Portanto, após V-1 rodadas, todos os caminhos devem estar ótimos.
    // 3. Detecção de Ciclo Negativo: Se na V-ésima rodada ainda for possível relaxar
    //    alguma aresta, isso prova matematicamente a existência de um ciclo negativo
    //    alcançável.
    // 4. Tolerância Numérica: Ao lidar com 'double', é vital usar um epsilon (1e-9)
    //    nas comparações para evitar falsos positivos de ciclos negativos.
    // ----------------------------------------------------------
    cout << "\n------- Bellman-Ford SP (s=" << s << ") -------" << endl;
    Process::BellmanFordSP bfSP(G, s);
    
    // Verificação de Ciclo Negativo (Resultado do algoritmo)
    if (bfSP.hasNegativeCycle()) {
        cout << "Ciclo negativo detectado! Arestas do ciclo:" << endl;
        stack<Edge> cycle = bfSP.getNegativeCycle();
        double cycleWeight = 0;
        while (!cycle.empty()) {
            Edge e = cycle.top(); cycle.pop();
            cout << e.from() << " -> " << e.to() << " (" << e.weight << ")" << endl;
            cycleWeight += e.weight;
        }
        cout << "Peso total do ciclo: " << cycleWeight << endl;
    } 
    else {
        for (int v = 0; v < G.getV(); v++) {
            if (bfSP.hasPathTo(v))
                 cout << "Caminho de " << s << " para " << v << " (" << fixed << setprecision(2) << bfSP.distanceTo(v) << ")" << endl;
        }    
    }

    // ==========================================================
    // 4. FLOYD-WARSHALL (Todos para Todos)
    // ----------------------------------------------------------
    // TIPO: Programação Dinâmica.
    // USO: Encontrar a menor distância entre TODOS os pares de vértices (All-Pairs).
    // COMPLEXIDADE: O(V^3) - Cúbica (Lento, usar apenas para grafos pequenos/densos).
    //
    // PECULIARIDADES:
    // 1. Matriz de Adjacência: Usa uma matriz D[i][j] em vez de listas.
    // 2. Lógica dos 3 Loops:
    //    - k: vértice intermediário permitido.
    //    - i: origem.
    //    - j: destino.
    //    - D[i][j] = min(D[i][j], D[i][k] + D[k][j]).
    // 3. Ciclos Negativos: Detectados checando a diagonal principal. Se dist[i][i] < 0,
    //    então o vértice 'i' faz parte de um ciclo negativo.
    // ----------------------------------------------------------
    cout << "\n------- Floyd-Warshall -------" << endl;
    Process::FloydWarshall fw(G);
    if(fw.hasNegCycle()) {
        cout << "Ciclo negativo detectado!" << endl;
    } else {
        cout << "      ";
        for(int i=0; i<G.getV(); i++) cout << setw(6) << i;
        cout << endl;
        for (int i = 0; i < G.getV(); i++) {
            cout << setw(3) << i << ": ";
            for (int j = 0; j < G.getV(); j++) {
                double d = fw.distTo(i, j);
                if (d > 99999) cout << setw(6) << "Inf";
                else cout << setw(6) << fixed << setprecision(2) << d;
            }
            cout << endl;
        }
    }

    // ---------------------------------------------------
    // 5. Execução do FORD-FULKERSON (Fluxo Máximo)
    // ---------------------------------------------------
    // TIPO: Método de Caminhos Aumentantes (Implementação Edmonds-Karp via BFS).
    // USO: Encontrar o fluxo máximo de uma fonte 's' para um sorvedouro 't'.
    // COMPLEXIDADE: O(V * E^2).
    //
    // PECULIARIDADES:
    // 1. Grafo Residual: O algoritmo precisa navegar "para trás" nas arestas para desfazer
    //    fluxos se necessário. Por isso, criamos arestas virtuais reversas.
    // 2. Capacidades: Assume capacidades não-negativas nas arestas.
    // 3. Gargalo (Bottleneck): O fluxo adicionado em cada iteração é limitado pela
    //    aresta de menor capacidade residual no caminho encontrado.
    // 4. Min-Cut (Dualidade): Ao terminar, os vértices que ainda são alcançáveis a partir
    //    de 's' no grafo residual formam o conjunto do Corte Mínimo.
    // ---------------------------------------------------
    int t = G.getV() - 1; // Define último vértice como destino (padrão para testes)
    cout << "\n------- Max Flow (Ford-Fulkerson s=" << s << " t=" << t << ") -------" << endl;
    Process::FordFulkerson maxflow(G, s, t);
    
    cout << "Fluxo Maximo: " << maxflow.getValue() << endl;
    cout << "Fluxo nas arestas:" << endl;
    for(int v=0; v < G.getV(); v++) {
        for(const auto& e : G.getAdj(v)) {
            // Exibe apenas arestas originais que estão sendo usadas pelo fluxo
            if(e.flow > 0 && e.capacity() > 0) 
                 cout << "  " << e.from() << " -> " << e.to() << " : " << e.flow << "/" << e.capacity() << endl;
        }
    }
    cout << "Vertices no Min-Cut: { ";
    // O corte mínimo é formado pelos vértices marcados (alcançáveis) na última BFS no grafo residual
    for(int v=0; v < G.getV(); v++) {
        if(maxflow.inCut(v)) cout << v << " ";
    }
    cout << "}" << endl;

    return 0;
}