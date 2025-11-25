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
    // 1. DIJKSTRA
    // ==========================================================
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
                path.push(e);
            return path;
        }       
    };

    // ==========================================================
    // 2. ACYCLIC SP
    // ==========================================================
    class AcyclicSP
    {
    private:
        vector<double> distTo;
        vector<Edge> edgeTo;
        bool cycleFound; 

    public:
        AcyclicSP(const EWDigraph &G, int s) : cycleFound(false)
        {
            distTo.resize(G.getV(), numeric_limits<double>::infinity());
            edgeTo.resize(G.getV(), Edge(-1, -1, 0.0));
            distTo[s] = 0.0;

            vector<bool> marked(G.getV(), false);
            vector<bool> onStack(G.getV(), false);
            stack<int> reversePost;
            
            function<void(int)> dfs = [&](int v) {
                marked[v] = true;
                onStack[v] = true;
                for (const Edge &e : G.getAdj(v)) {
                    if (cycleFound) return;
                    int w = e.to();
                    if (!marked[w]) dfs(w);
                    else if (onStack[w]) cycleFound = true;
                }
                onStack[v] = false;
                reversePost.push(v);
            };

            for (int v = 0; v < G.getV(); v++)
                if (!marked[v] && !cycleFound) dfs(v);

            if (cycleFound) return;

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
    // 3. BELLMAN-FORD SP
    // ==========================================================
    class BellmanFordSP
    {
    private:
        vector<double> distTo;
        vector<Edge> edgeTo;
        stack<Edge> cycle;
        bool hasCycle;

    public:
        BellmanFordSP(const EWDigraph &G, int s) : hasCycle(false)
        {
            distTo.resize(G.getV(), numeric_limits<double>::infinity());
            edgeTo.resize(G.getV(), Edge(-1, -1, 0.0));
            distTo[s] = 0.0;    

            for (int i = 0; i < G.getV() - 1; i++) {
                for (int v = 0; v < G.getV(); v++) {
                    for (const Edge &e : G.getAdj(v)) {
                        relax(e);
                    }
                }
            }

            for (int v = 0; v < G.getV(); v++) {
                for (const Edge &e : G.getAdj(v)) {
                    int w = e.to();
                    // Tolerância de 1e-9 para evitar falso positivo com double
                    if (distTo[v] != numeric_limits<double>::infinity() && 
                        distTo[w] > distTo[v] + e.weight + 1e-9)
                    {
                        hasCycle = true;
                        
                        int cycleNode = w;
                        for(int i = 0; i < G.getV(); ++i) {
                           if (edgeTo[cycleNode].from() == -1) break; 
                           cycleNode = edgeTo[cycleNode].from();
                        }

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
        stack<Edge> getNegativeCycle() const { return cycle; }
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
    // 4. FLOYD-WARSHALL (Nova Implementação)
    // ==========================================================
    class FloydWarshall
    {
    private:
        vector<vector<double>> dist;
        bool hasNegativeCycle;

    public:
        FloydWarshall(const EWDigraph &G) : hasNegativeCycle(false)
        {
            int V = G.getV();
            dist.resize(V, vector<double>(V, numeric_limits<double>::infinity()));

            for (int v = 0; v < V; v++) {
                dist[v][v] = 0.0;
                for (const Edge &e : G.getAdj(v)) {
                    dist[e.from()][e.to()] = e.weight;
                }
            }

            for (int k = 0; k < V; k++) {
                for (int i = 0; i < V; i++) {
                    for (int j = 0; j < V; j++) {
                        if (dist[i][k] != numeric_limits<double>::infinity() &&
                            dist[k][j] != numeric_limits<double>::infinity()) {
                            if (dist[i][j] > dist[i][k] + dist[k][j] + 1e-9) {
                                dist[i][j] = dist[i][k] + dist[k][j];
                            }
                        }
                    }
                }
                if (dist[k][k] < -1e-9) { hasNegativeCycle = true; return; }
            }
        }

        double distTo(int s, int t) { return dist[s][t]; }
        bool hasNegCycle() { return hasNegativeCycle; }
    };

    // ==========================================================
    // 5. FORD-FULKERSON (Nova Implementação)
    // ==========================================================
    class FordFulkerson
    {
    private:
        vector<bool> marked;
        vector<Edge*> edgeTo;
        double value;

        bool hasAugmentingPath(EWDigraph &G, int s, int t) {
            edgeTo.assign(G.getV(), nullptr);
            marked.assign(G.getV(), false);
            queue<int> q;
            q.push(s);
            marked[s] = true;

            while (!q.empty() && !marked[t]) {
                int v = q.front(); q.pop();
                for (Edge &e : G.getAdj(v)) {
                    int w = e.other(v);
                    if (e.residualCapacityTo(w) > 0 && !marked[w]) {
                        edgeTo[w] = &e;
                        marked[w] = true;
                        q.push(w);
                    }
                }
            }
            return marked[t];
        }

    public:
        FordFulkerson(EWDigraph &G, int s, int t) : value(0.0) {
            // Adiciona arestas reversas necessárias para o grafo residual
            vector<Edge> originalEdges = G.getAllEdges();
            for(const auto& e : originalEdges) {
                bool hasReverse = false;
                for(const auto& rev : G.getAdj(e.to())) {
                    if(rev.to() == e.from()) { hasReverse = true; break; }
                }
                if(!hasReverse) G.addEdge(Edge(e.to(), e.from(), 0.0));
            }

            while (hasAugmentingPath(G, s, t)) {
                double bottle = numeric_limits<double>::infinity();
                for (int v = t; v != s; v = edgeTo[v]->other(v)) {
                    bottle = min(bottle, edgeTo[v]->residualCapacityTo(v));
                }
                for (int v = t; v != s; v = edgeTo[v]->other(v)) {
                    edgeTo[v]->addResidualFlowTo(v, bottle);
                }
                value += bottle;
            }
        }
        double getValue() const { return value; }
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
    G.showDot();
    
    int s = 0; 

    // 1. Dijkstra
    cout << "\n------- Dijkstra SP (s=" << s << ") -------" << endl;
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
                    cout << e.from() << "->" << e.to() << " ";
                }
                cout << endl;
            } else {
                cout << "Nao ha caminho de " << s << " para " << v << endl;
            }
        }
    }

    // 2. Acyclic SP
    cout << "\n------- Acyclic SP (s=" << s << ") -------" << endl;
    Process::AcyclicSP acyclicSP(G, s);
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

    // 3. Bellman-Ford
    cout << "\n------- Bellman-Ford SP (s=" << s << ") -------" << endl;
    Process::BellmanFordSP bfSP(G, s);
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

    // 4. Floyd-Warshall
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

    // 5. Ford-Fulkerson
    int t = G.getV() - 1; // Define último vértice como destino
    cout << "\n------- Max Flow (Ford-Fulkerson s=" << s << " t=" << t << ") -------" << endl;
    Process::FordFulkerson maxflow(G, s, t);
    
    cout << "Fluxo Maximo: " << maxflow.getValue() << endl;
    cout << "Fluxo nas arestas:" << endl;
    for(int v=0; v < G.getV(); v++) {
        for(const auto& e : G.getAdj(v)) {
            if(e.flow > 0 && e.capacity() > 0) 
                 cout << "  " << e.from() << " -> " << e.to() << " : " << e.flow << "/" << e.capacity() << endl;
        }
    }
    cout << "Vertices no Min-Cut: { ";
    for(int v=0; v < G.getV(); v++) {
        if(maxflow.inCut(v)) cout << v << " ";
    }
    cout << "}" << endl;

    return 0;
}