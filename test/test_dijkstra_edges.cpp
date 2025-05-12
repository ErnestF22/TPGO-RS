#include <iostream>
#include <vector>
#include <queue>
#include <set>
#include <limits.h>

std::vector<bool> findAnswer(int n, std::vector<std::vector<int>> &edges)
{
    std::vector<std::pair<int, int>> adj[n];
    int totalEdges = edges.size();
    for (int i = 0; i < totalEdges; i++)
    {
        int u = edges[i][0];
        int v = edges[i][1];
        int w = edges[i][2];
        adj[u].push_back({v, w});
        adj[v].push_back({u, w});
    }

    std::vector<int> dist(n, INT_MAX);
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>> pq;
    pq.push({0, 0});
    dist[0] = 0;

    while (!pq.empty())
    {
        int node = pq.top().second;
        int d = pq.top().first;
        pq.pop();

        for (auto it : adj[node])
        {
            int v = it.first;
            int w = it.second;
            if (dist[v] > dist[node] + w)
            {
                dist[v] = dist[node] + w;
                pq.push({dist[v], v});
            }
        }
    }

    std::vector<bool> ans(totalEdges, false);

    if (dist[n - 1] == INT_MAX)
        return ans;

    std::vector<int> vis(n, 0);
    std::queue<std::pair<int, int>> q;
    q.push({n - 1, dist[n - 1]});

    std::set<std::pair<int, int>> st;
    while (!q.empty())
    {
        int node = q.front().first;
        int d = q.front().second;
        q.pop();
        vis[node] = 1;
        for (auto it : adj[node])
        {
            int v = it.first;
            int w = it.second;
            int distRemain = dist[node] - w;
            if (distRemain == dist[v])
            {
                st.insert({node, v});
                st.insert({v, node});
                q.push({v, dist[v]});
            }
        }
    }

    for (int i = 0; i < totalEdges; i++)
    {
        int u = edges[i][0];
        int v = edges[i][1];
        if (st.find({u, v}) != st.end())
            ans[i] = true;
    }

    return ans;
}

int main(int argc, char **argv)
{

    std::vector<std::vector<int>> edges;
    edges.push_back(std::vector<int>({1, 0, 1}));
    edges.push_back(std::vector<int>({2, 0, 1}));
    edges.push_back(std::vector<int>({3, 0, 1}));
    edges.push_back(std::vector<int>({0, 1, 1}));
    edges.push_back(std::vector<int>({2, 1, 1}));
    edges.push_back(std::vector<int>({3, 1, 1}));
    edges.push_back(std::vector<int>({4, 1, 1}));
    edges.push_back(std::vector<int>({0, 2, 1}));
    edges.push_back(std::vector<int>({1, 2, 1}));
    edges.push_back(std::vector<int>({3, 2, 1}));
    edges.push_back(std::vector<int>({4, 2, 1}));
    edges.push_back(std::vector<int>({0, 3, 1}));
    edges.push_back(std::vector<int>({1, 3, 1}));
    edges.push_back(std::vector<int>({2, 3, 1}));
    edges.push_back(std::vector<int>({4, 3, 1}));
    edges.push_back(std::vector<int>({1, 4, 1}));
    edges.push_back(std::vector<int>({2, 4, 1}));
    edges.push_back(std::vector<int>({3, 4, 1}));

    auto asd = findAnswer(5, edges);

    for (int i = 0; i < asd.size(); ++i)
        std::cout << asd[i] << std::endl;

    return 0;
}