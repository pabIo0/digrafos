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

class AcyclicSP
{
private:
    vector<double> distTo;
    vector<Edge> edgeTo;
public:
    AcyclicSP(const EWDigraph &G, int s)
    {
        distTo.resize(G.getV(), numeric_limits<double>::infinity());
        edgeTo.resize(G.getV(), Edge(-1, -1, 0.0));
        distTo[s] = 0.0;    
        // Ordenação topológica
        vector<bool> marked(G.getV(), false);
        stack<int> reversePost;
        function<void(int)> dfs = [&](int v) {
            marked[v] = true;
            for (const Edge &e : G.getAdj(v))
            {
                int w = e.to();
                if (!marked[w]) dfs(w);
            }
            reversePost.push(v);
        };
        for (int v = 0; v < G.getV(); v++)
        {
            if (!marked[v]) dfs(v);
        }
        // Relaxamento das arestas na ordem topológica      
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
public:
    belmanFordSP(const EWDigraph &G, int s)
    {
        distTo.resize(G.getV(), numeric_limits<double>::infinity());
        edgeTo.resize(G.getV(), Edge(-1, -1, 0.0));
        distTo[s] = 0.0;    
        // Relaxamento das arestas |V|-1 vezes
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

    Process::AcyclicSP acyclicSP(G, 0); // Caminhos mínimos em grafo acíclico a partir do vértice 0
    cout << " ------- Acyclic SP -------" << endl;
    for (int v = 0; v < G.getV(); v++)      
    {
        if (acyclicSP.hasPathTo(v))
        {
            stack<Edge> path = acyclicSP.pathTo(v);
            cout << "Caminho minimo de 0 para " << v << " (distancia: " << fixed << setprecision(2) << acyclicSP.distanceTo(v) << "): ";
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

    Process::belmanFordSP bfSP(G, 0); // Caminhos mínimos usando Bellman-Ford a partir do vértice 0
    cout << " ------- Bellman-Ford SP -------" << endl;
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
    cout << endl;


    return 0;
}

