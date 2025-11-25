#include "EWDigraph.h"
#include <functional>

using namespace std;

class Process
{
    private:
        EWDigraph G;

    public:
        Process(const EWDigraph &G) : G(G) {}

class DjkstraSP
{
    private:
        vector<double> distTo;
        vector<Edge> edgeTo;
        IndexMinPQ<double> pq;      
    public:
        DjkstraSP(const EWDigraph &G, int s) : distTo(G.getV(), numeric_limits<double>::infinity()), edgeTo(G.getV()), pq(G.getV())
        {
            distTo[s] = 0.0;    
            pq.insert(s, 0.0);
            while (!pq.isEmpty())
            {
                int v = pq.delMin();
                for (const Edge &e : G.getAdj(v))
                {
                    relax(e);
                }
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
                if (pq.contains(w)) pq.decreaseKey(w, distTo[w]);
                else pq.insert(w, distTo[w]);
            }
        }
        double distanceTo(int v) const { return distTo[v]; }
        bool hasPathTo(int v) const { return distTo[v] < numeric_limits<double>::infinity(); }
        stack<Edge> pathTo(int v) const
        {   
            stack<Edge> path;
            if (!hasPathTo(v)) return path;
            for (Edge e = edgeTo[v]; e.from() != -1; e = edgeTo[e.from()])
            {
                path.push(e);
            }
            return path;
        }       
    };

    private:
        bool isCyclicAux(int v, vector<bool>& visited, vector<bool>& onStack) {
            visited[v] = true;
            onStack[v] = true;

            for (const Edge &e : G.getAdj(v)) {
                int w = e.to();
                if (!visited[w]) {
                    if (isCyclicAux(w, visited, onStack)) return true;
                } else if (onStack[w]) {
                    return true; // Ciclo detectado!
                }
            }
            onStack[v] = false;
            return false;
    }

    public:
        // Função pública que verifica se o grafo tem ciclo
        bool hasCycle() {
            vector<bool> visited(G.getV(), false);
            vector<bool> onStack(G.getV(), false);
            for (int v = 0; v < G.getV(); v++) {
                if (!visited[v]) {
                    if (isCyclicAux(v, visited, onStack)) return true;
                }
            }
            return false;
    }

class AcyclicSP
    {
    private:
        vector<double> distTo;
        vector<Edge> edgeTo;
        bool cycleFound; // Flag para indicar se ciclo foi encontrado

    public:
        AcyclicSP(const EWDigraph &G, int s) : cycleFound(false) // Inicializa flag
        {
            distTo.resize(G.getV(), numeric_limits<double>::infinity());
            edgeTo.resize(G.getV(), Edge(-1, -1, 0.0));
            distTo[s] = 0.0;

            // Estruturas para Ordenação Topológica e Detecção de Ciclo
            vector<bool> marked(G.getV(), false);
            vector<bool> onStack(G.getV(), false); // Rastreia a pilha de recursão atual
            stack<int> reversePost;
            
            // Função Lambda DFS com Detecção de Ciclo
            function<void(int)> dfs = [&](int v) {
                marked[v] = true;
                onStack[v] = true; // Vértice entra na pilha de recursão

                for (const Edge &e : G.getAdj(v))
                {
                    int w = e.to();
                    if (cycleFound) return; // Se já achou ciclo, aborta
                    
                    if (!marked[w]) {
                        dfs(w);
                    } else if (onStack[w]) {
                        // Se o vizinho já foi visitado E está na pilha atual -> CICLO!
                        cycleFound = true;
                    }
                }
                
                onStack[v] = false; // Vértice sai da pilha de recursão
                reversePost.push(v);
            };

            // Executa a DFS para todos os componentes
            for (int v = 0; v < G.getV(); v++)
            {
                if (!marked[v] && !cycleFound) dfs(v);
            }

            // Se encontrou ciclo, não faz relaxamento e avisa (ou lança exceção)
            if (cycleFound) {
                // Apenas retorna, deixando distTo como infinito (exceto a origem)
                // O main deve checar hasCycle() antes de confiar nos resultados
                return; 
            }

            // Se não tem ciclo, prossegue com o relaxamento na ordem topológica
            while (!reversePost.empty())
            {
                int v = reversePost.top();
                reversePost.pop();
                for (const Edge &e : G.getAdj(v))
                {
                    relax(e);
                }
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

        // Método público para verificar se ciclo foi detectado
        bool hasCycle() const { return cycleFound; }
        bool isValid() const { return !cycleFound; }
        double distanceTo(int v) const { return distTo[v]; }
        bool hasPathTo(int v) const { return distTo[v] < numeric_limits<double>::infinity(); }
        stack<Edge> pathTo(int v) const
        {
            stack<Edge> path;
            if (!hasPathTo(v)) return path;
            for (Edge e = edgeTo[v]; e.from() != -1; e = edgeTo[e.from()])
            {
                path.push(e);
            }
            return path;
        }
    };

class belmanFordSP
    {
    private:
        vector<double> distTo;
        vector<Edge> edgeTo;
        bool hasNegCycle; // 1. Nova variável para armazenar o estado

    public:
        belmanFordSP(const EWDigraph &G, int s) : hasNegCycle(false) // Inicializa como false
        {
            distTo.resize(G.getV(), numeric_limits<double>::infinity());
            edgeTo.resize(G.getV(), Edge(-1, -1, 0.0));
            distTo[s] = 0.0;    

            // Passo 1: Relaxamento das arestas |V|-1 vezes (Seu código original)
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

            // Passo 2: Verificar Ciclos Negativos
            // Passamos por todas as arestas mais uma vez.
            for (int v = 0; v < G.getV(); v++)
            {
                for (const Edge &e : G.getAdj(v))
                {
                    int w = e.to();
                    // Se ainda conseguirmos melhorar a distância, existe um ciclo negativo
                    if (distTo[v] != numeric_limits<double>::infinity() && 
                        distTo[w] > distTo[v] + e.weight)
                    {
                        hasNegCycle = true;
                        return; // Encontramos, podemos parar o construtor
                    }
                }
            }
        }   

        void relax(const Edge &e)
        {
            int v = e.from();
            int w = e.to();
            // Verificação extra para evitar somar com infinito
            if (distTo[v] != numeric_limits<double>::infinity() && 
                distTo[w] > distTo[v] + e.weight)
            {
                distTo[w] = distTo[v] + e.weight;
                edgeTo[w] = e;
            }
        }

        // Métodos auxiliares
        bool hasNegativeCycle() const { return hasNegCycle; } // Método para o main consultar
        double distanceTo(int v) const { return distTo[v]; }
        bool hasPathTo(int v) const { return distTo[v] < numeric_limits<double>::infinity(); }
        
        stack<Edge> pathTo(int v) const
        {   
            stack<Edge> path;
            if (!hasPathTo(v)) return path;
            for (Edge e = edgeTo[v]; e.from() != -1; e = edgeTo[e.from()])
            {
                path.push(e);
            }
            return path;
        }
    };
};

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        cerr << "Uso: " << argv[0] << " <arquivo_do_grafo>" << endl;
        return 1;
    }

    ifstream in(argv[1]);
    if (!in)
    {
        cerr << "Erro ao abrir o arquivo: " << argv[1] << endl;
        return 1;
    }

    // Cria o grafo a partir do arquivo
    EWDigraph G(in);

    G.show();
    G.showDot();
    
    Process::DjkstraSP dsp(G, 0); // Caminhos mínimos a partir do vértice 0
    cout << " ------- Dijkstra -------" << endl;
    for (int v = 0; v < G.getV(); v++)
    {
        if (dsp.hasPathTo(v))
        {
            stack<Edge> path = dsp.pathTo(v);
            cout << "Caminho minimo de 0 para " << v << " (distancia: " << fixed << setprecision(2) << dsp.distanceTo(v) << "): ";
            while (!path.empty())
            {
                Edge e = path.top();
                cout << e.from() << "->" << e.to() << " (" << fixed << setprecision(2) << e.weight << ") ";
                path.pop();
            }
            cout << endl;
        }
        else
        {
            cout << "Nao ha caminho de 0 para " << v << endl;
        }
    }
    cout << endl;

    Process::AcyclicSP acyclicSP(G, 0); 
    cout << " ------- Acyclic SP -------" << endl;
    
    if (acyclicSP.isValid()) {
        for (int v = 0; v < G.getV(); v++)      
        {
            if (acyclicSP.hasPathTo(v)) {
                stack<Edge> path = acyclicSP.pathTo(v);
                cout << "Caminho minimo de 0 para " << v << " (distancia: " << fixed << setprecision(2) << acyclicSP.distanceTo(v) << "): ";
                while (!path.empty())
                {
                    Edge e = path.top();
                    cout << e.from() << "->" << e.to() << " (" << fixed << setprecision(2) << e.weight << ") ";
                    path.pop();
                }
                cout << endl;
            } else {
                cout << "Nao ha caminho de 0 para " << v << endl;
            }
        }
    } else {
        // --- ADICIONE ESTE BLOCO ELSE ---
        cout << "Erro: O grafo contem ciclos! O algoritmo de caminho minimo aciclico nao pode ser aplicado." << endl;
    }
    cout << endl;

    Process::belmanFordSP bfSP(G, 0); 
    cout << " ------- Bellman-Ford SP -------" << endl;
    
    // Verifica primeiro se existe ciclo negativo
    if (bfSP.hasNegativeCycle()) {
        cout << "Ciclo negativo detectado! Nao e possivel calcular caminhos minimos confiaveis." << endl;
    } 
    else {
        for (int v = 0; v < G.getV(); v++)      
        {
            if (bfSP.hasPathTo(v))
            {
                stack<Edge> path = bfSP.pathTo(v);
                cout << "Caminho minimo de 0 para " << v << " (distancia: " << fixed << setprecision(2) << bfSP.distanceTo(v) << "): ";
                while (!path.empty())
                {
                    Edge e = path.top();
                    cout << e.from() << "->" << e.to() << " (" << fixed << setprecision(2) << e.weight << ") ";
                    path.pop();
                }
                cout << endl;
            }
            else
            {
                cout << "Nao ha caminho de 0 para " << v << endl;
            }
        }   
    }
    cout << endl;

    return 0;
}