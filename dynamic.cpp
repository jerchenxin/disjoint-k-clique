#include <omp.h>

#include "dynamic.h"

DynamicGraph::DynamicGraph(string fileName, int k) : k(k) {
    BuildIndex(fileName);
}

void DynamicGraph::BuildIndex(string fileName) {
    ThirdGreedy* thirdGreedy = new ThirdGreedy(fileName, k);
    thirdGreedy->solve(fileName);
    output = move(thirdGreedy->output);
    graph = move(thirdGreedy->graph);
    vertexRank = move(thirdGreedy->vertexRank);
    n = thirdGreedy->n;
    reverseGraph.resize(n);

    for (int u=0;u<n;u++) {
        for (auto v : graph[u]) {
            reverseGraph[v].emplace_back(u);
        }
    }

    delete thirdGreedy;

    auto start = std::chrono::high_resolution_clock::now();
    BuildDynamicIndex();
    auto end = std::chrono::high_resolution_clock::now();
    auto diff = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    cout << endl <<  "Index Time: " << diff.count() << endl;

    Info();
}

void DynamicGraph::BuildDynamicIndex() {
    group.resize(n, -1);
    index.resize(n);

    for (auto& c : output) {
        for (auto i : c) {
            group[i] = c.front();
        }

        invGroup[c.front()] = c;
    }

    subgraph = KSubgraph(k, n);

    // for (auto& c : output) {
    //     FindIndexByClique(c);
    // }

    omp_set_num_threads(64);
    int i;
    int cliqueNum = output.size();
    #pragma omp parallel private(i)
    {
        vector<set<Clique*>> tmpIndex(n);
        int pid = omp_get_thread_num();
        KSubgraph subgraph(k, n);
        #pragma omp for schedule(dynamic, 1) nowait
        for (i=0;i<cliqueNum;i++) {
            MakeSubGraph(subgraph, k, output[i]);

            for (auto i : subgraph.nodes[k]) {
                if (subgraph.degree[k][i] < k-1) {
                    continue;
                }

                BuildSubGraph(subgraph, k, i);
                vector<int> prev = {i};
                FindAllCliqueMulti(subgraph, k-1, prev, tmpIndex);
            }
        }

        #pragma omp critical
        {
            for (int j=0;j<n;j++) {
                index[j].insert(tmpIndex[j].begin(), tmpIndex[j].end());
            }
        }
    }
}

void DynamicGraph::FindIndexByClique(vector<int>& c) {
    MakeSubGraph(subgraph, k, c);

    for (auto i : subgraph.nodes[k]) {
        if (subgraph.degree[k][i] < k-1) {
            continue;
        }

        BuildSubGraph(subgraph, k, i);
        vector<int> prev = {i};
        FindAllClique(subgraph, k-1, prev);
    }
}

void DynamicGraph::MakeSubGraph(KSubgraph& subgraph, int k, vector<int>& c) {
    subgraph.nodes[k].clear();

    // mark
    for (auto i : c) {
        subgraph.nodes[k].emplace(i);
        subgraph.candidate[i] = 1;

        for (auto v : graph[i]) {
            if (group[v] == -1 || group[v] == group[i]) {
                subgraph.nodes[k].emplace(v);
                subgraph.candidate[v] = 1;
            }
        }

        for (auto v : reverseGraph[i]) {
            if (group[v] == -1 || group[v] == group[i]) {
                subgraph.nodes[k].emplace(v);
                subgraph.candidate[v] = 1;
            }
        }
    }

    // extract edge
    for (auto i : subgraph.nodes[k]) {
        subgraph.graph[i].clear();
        for (auto v : graph[i]) {
            if (subgraph.candidate[v] == 1) {
                subgraph.graph[i].emplace_back(v);
            }
        }
    }

    for (auto i : subgraph.nodes[k]) {
        subgraph.degree[k][i] = subgraph.graph[i].size();
    }

    for (auto i : subgraph.nodes[k]) {
        subgraph.candidate[i] = 0;
    }
}

void DynamicGraph::BuildSubGraph(KSubgraph& subgraph, int k, int u) {
    auto& nodes = subgraph.nodes;
    auto& degree = subgraph.degree;
    auto& candidate = subgraph.candidate;
    auto& graph = subgraph.graph;

    nodes[k-1].clear();

    for (int i=0;i<degree[k][u];i++) {
        int v = graph[u][i];
        candidate[v] = 1;
        nodes[k-1].emplace(v);
    }

    for (int i=0;i<degree[k][u];i++) {
        int v = graph[u][i];
        int nextIndex = 0;
        for (int j=0;j<degree[k][v];j++) {
            if (candidate[graph[v][j]] == 1) {
                swap(graph[v][j], graph[v][nextIndex]);
                nextIndex++;
            }
        }
        degree[k-1][v] = nextIndex;
    }

    for (int i=0;i<degree[k][u];i++) {
        int v = graph[u][i];
        candidate[v] = 0;
    }
}

void DynamicGraph::FindAllClique(KSubgraph& subgraph, int k, vector<int>& prev) {
    auto& nodes = subgraph.nodes;
    auto& degree = subgraph.degree;
    auto& graph = subgraph.graph;

    if (k == 2) {
        for (auto i : nodes[k]) {
            for (int j=0;j<degree[k][i];j++) {
                vector<int> tmp = prev;
                tmp.push_back(i);
                tmp.push_back(graph[i][j]);

                for (auto p : tmp) {
                    if (group[p] == -1) {
                        Clique* c = new Clique();
                        c->nodes = move(tmp);

                        for (auto s : c->nodes) {
                            index[s].emplace(c);
                        }
                        break;
                    }
                }
            }
        }
    } else {
        for (auto i : nodes[k]) {
            if (degree[k][i] < k-1) {
                continue;
            }

            BuildSubGraph(subgraph, k, i);
            prev.push_back(i);
            FindAllClique(subgraph, k-1, prev);
            prev.pop_back();
        }
    }
}

void DynamicGraph::FindAllCliqueMulti(KSubgraph& subgraph, int k, vector<int>& prev, vector<set<Clique*>>& tmpIndex) {
    auto& nodes = subgraph.nodes;
    auto& degree = subgraph.degree;
    auto& graph = subgraph.graph;

    if (k == 2) {
        for (auto i : nodes[k]) {
            for (int j=0;j<degree[k][i];j++) {
                vector<int> tmp = prev;
                tmp.push_back(i);
                tmp.push_back(graph[i][j]);

                for (auto p : tmp) {
                    if (group[p] == -1) {
                        Clique* c = new Clique();
                        c->nodes = move(tmp);

                        for (auto s : c->nodes) {
                            tmpIndex[s].emplace(c);
                        }
                        break;
                    }
                }
            }
        }
    } else {
        for (auto i : nodes[k]) {
            if (degree[k][i] < k-1) {
                continue;
            }

            BuildSubGraph(subgraph, k, i);
            prev.push_back(i);
            FindAllCliqueMulti(subgraph, k-1, prev, tmpIndex);
            prev.pop_back();
        }
    }
}

// do not allow multi edge
void DynamicGraph::Insertion(int u, int v) {
    AddEdge(u, v);

    // group node -- group node: pass
    if (group[u] != group[v] && group[v] != -1 && group[u] != -1) {
        return;
    }
    
    if (group[u] != group[v]) { // free node -- group node
        // update index
        if (group[u] == -1) {
            swap(u, v);
        }

        set<int> modifiedGroup;

        MakeSubGraph(subgraph, k, u, v);

        for (auto i : subgraph.nodes[k]) {
            if (subgraph.degree[k][i] < k-1) {
                continue;
            }

            BuildSubGraph(subgraph, k, i);
            vector<int> prev = {i};
            FindAllClique(subgraph, k-1, prev, modifiedGroup, u, v);
        }

        // index changes, find the modified group
        Swap(modifiedGroup);
    } else { // free node -- free node
        // find free clique
        set<Clique*> freeCliqueList;

        MakeFreeSubgraph(subgraph, k, u, v);

        for (auto i : subgraph.nodes[k]) {
            if (subgraph.degree[k][i] < k-1) {
                continue;
            }

            BuildSubGraph(subgraph, k, i);
            vector<int> prev = {i};
            FindAllFreeClique(subgraph, k-1, prev, freeCliqueList);
        }

        auto result = DisjointClique(freeCliqueList); 
        if (!result.empty()) {
            set<int> newGroup;

            // update group
            for (auto targetC : result) {
                for (auto node : targetC) {
                    group[node] = targetC.front();
                }
                invGroup[targetC.front()] = targetC;
                newGroup.emplace(targetC.front());
            }

            // rm invalid index
            for (auto id : newGroup) {
                for (auto node : invGroup[id]) {
                    auto nodeIndex = move(index[node]);
                    for (auto i : nodeIndex) {
                        DeleteIndex(i);
                    }
                }
            }

            // update index
            for (auto id : newGroup) {
                FindIndexByClique(invGroup[id]);
            }
        } else {
            // update index
            set<int> modifiedGroup;
            MakeSubgraphByFreeEdge(subgraph, k, u, v);

            for (auto i : subgraph.nodes[k]) {
                if (subgraph.degree[k][i] < k-1) {
                    continue;
                }

                BuildSubGraph(subgraph, k, i);
                vector<int> prev = {i};
                FindAllCliqueByFreeEdge(subgraph, k-1, prev, modifiedGroup, u, v);
            }

            // try swap
            Swap(modifiedGroup);
        }
    }
}

void DynamicGraph::AddEdge(int u, int v) {
    if (vertexRank[u] < vertexRank[v]) {
        swap(u, v);
    }

    graph[u].push_back(v);
    reverseGraph[v].push_back(u);
}

// s should be the group node
void DynamicGraph::MakeSubGraph(KSubgraph& subgraph, int k, int s, int t) {
    subgraph.nodes[k].clear();

    // mark
    subgraph.nodes[k].emplace(s);
    subgraph.candidate[s] = 1;
    subgraph.nodes[k].emplace(t);
    subgraph.candidate[t] = 1;

    set<int> candS;

    for (auto v : graph[s]) {
        if (group[v] == -1 || group[v] == group[s]) {
            candS.emplace(v);
        }
    }

    for (auto v : reverseGraph[s]) {
        if (group[v] == -1 || group[v] == group[s]) {
            candS.emplace(v);
        }
    }

    for (auto v : graph[t]) {
        if (candS.find(v) != candS.end()) {
            subgraph.nodes[k].emplace(v);
            subgraph.candidate[v] = 1;
        }
    }

    for (auto v : reverseGraph[t]) {
        if (candS.find(v) != candS.end()) {
            subgraph.nodes[k].emplace(v);
            subgraph.candidate[v] = 1;
        }
    }

    // extract edge
    for (auto i : subgraph.nodes[k]) {
        subgraph.graph[i].clear();
        for (auto v : graph[i]) {
            if (subgraph.candidate[v] == 1) {
                subgraph.graph[i].emplace_back(v);
            }
        }
    }

    for (auto i : subgraph.nodes[k]) {
        subgraph.degree[k][i] = subgraph.graph[i].size();
    }

    for (auto i : subgraph.nodes[k]) {
        subgraph.candidate[i] = 0;
    }
}

void DynamicGraph::FindAllClique(KSubgraph& subgraph, int k, vector<int>& prev, set<int>& modifiedGroup, int s1, int t1) {
    auto& nodes = subgraph.nodes;
    auto& degree = subgraph.degree;
    auto& graph = subgraph.graph;

    if (k == 2) {
        for (auto i : nodes[k]) {
            for (int j=0;j<degree[k][i];j++) {
                vector<int> tmp = prev;
                tmp.push_back(i);
                tmp.push_back(graph[i][j]);

                set<int> nodeSet(tmp.begin(), tmp.end());

                // should contain the new edge
                if (nodeSet.find(s1)!=nodeSet.end() && nodeSet.find(t1)!=nodeSet.end()) {
                    Clique* c = new Clique();
                    c->nodes = move(tmp);

                    for (auto s : c->nodes) {
                        index[s].emplace(c);
                    }

                    for (auto s : c->nodes) {
                        if (group[s] != -1) {
                            modifiedGroup.emplace(group[s]);
                            break;
                        }
                    }
                }
            }
        }
    } else {
        for (auto i : nodes[k]) {
            if (degree[k][i] < k-1) {
                continue;
            }

            BuildSubGraph(subgraph, k, i);
            prev.push_back(i);
            FindAllClique(subgraph, k-1, prev, modifiedGroup, s1, t1);
            prev.pop_back();
        }
    }
}

 void DynamicGraph::Swap(set<int>& modifiedGroup) {
    queue<int> q;
    for (auto i : modifiedGroup) {
        q.push(i);
    }

    while (!q.empty()) {
        int groupID = q.front();
        q.pop();

        if (group[groupID] == -1) {
            continue;
        }

        modifiedGroup.clear();
        TrySwap(groupID, modifiedGroup);
        for (auto i : modifiedGroup) {
            q.push(i);
        }
    }
 }

void DynamicGraph::TrySwap(int groupID, set<int>& modifiedGroup) {
    set<Clique*> cs;
    for (auto i : invGroup[groupID]) {
        cs.insert(index[i].begin(), index[i].end());
    }

    auto result = DisjointClique(cs);
    if (result.size() > 1) {
        vector<int> nodes = move(invGroup[groupID]);
        invGroup.erase(groupID);

        // delete the index of this clique
        for (auto i : nodes) {
            group[i] = -1;

            auto nodeIndex = move(index[i]);
            for (auto i : nodeIndex) {
                DeleteIndex(i);
            }
        }

        set<int> newGroupID;

        for (auto targetC : result) {
            for (auto node : targetC) {
                group[node] = targetC.front();
            }
            invGroup[targetC.front()] = targetC;
            newGroupID.emplace(targetC.front());
        }

        // delete the index of the new cliques
        for (auto id : newGroupID) {
            for (auto node : invGroup[id]) {
                auto nodeIndex = move(index[node]);
                for (auto i : nodeIndex) {
                    DeleteIndex(i);
                }
            }
        }

        // update index of new cliques
        for (auto id : newGroupID) { 
            FindIndexByClique(invGroup[id]); // should modify
            modifiedGroup.emplace(id);
        }

        // update index of old cliques using the deleted clique
        set<int> newFreeNodes(nodes.begin(), nodes.end());
        for (auto targetC : result) {
            for (auto node : targetC) {
                newFreeNodes.erase(node);
            }
        }

        set<int> modifiedGroup(newGroupID.begin(), newGroupID.end());

        if (!newFreeNodes.empty()) {
            MakeSubgraphByFreeNodes(subgraph, k, newFreeNodes, modifiedGroup);

            for (auto i : subgraph.nodes[k]) {
                if (subgraph.degree[k][i] < k-1) {
                    continue;
                }

                BuildSubGraph(subgraph, k, i);
                vector<int> prev = {i};
                FindAllCliqueByFreeNodes(subgraph, k-1, prev, modifiedGroup, newFreeNodes);
            }
        }
    }
}

void DynamicGraph::MakeSubGraph(KSubgraph& subgraph, int k, int groupID) {
    subgraph.nodes[k].clear();

    // mark
    for (auto i : invGroup[groupID]) {
        subgraph.nodes[k].emplace(i);
        subgraph.candidate[i] = 1;
    }

    for (auto i : invGroup[groupID]) {
        for (auto v : graph[i]) {
            if (group[v] == -1 || group[v] == group[i]) {
                subgraph.nodes[k].emplace(v);
                subgraph.candidate[v] = 1;
            }
        }

        for (auto v : reverseGraph[i]) {
            if (group[v] == -1 || group[v] == group[i]) {
                subgraph.nodes[k].emplace(v);
                subgraph.candidate[v] = 1;
            }
        }
    }

    // extract edge
    for (auto i : subgraph.nodes[k]) {
        subgraph.graph[i].clear();
        for (auto v : graph[i]) {
            if (subgraph.candidate[v] == 1) {
                subgraph.graph[i].emplace_back(v);
            }
        }
    }

    for (auto i : subgraph.nodes[k]) {
        subgraph.degree[k][i] = subgraph.graph[i].size();
    }

    for (auto i : subgraph.nodes[k]) {
        subgraph.candidate[i] = 0;
    }
}

void DynamicGraph::MakeFreeSubgraph(KSubgraph& subgraph, int k, int s, int t) {
    subgraph.nodes[k].clear();

    // mark
    subgraph.nodes[k].emplace(s);
    subgraph.candidate[s] = 1;
    subgraph.nodes[k].emplace(t);
    subgraph.candidate[t] = 1;

    for (auto i : graph[s]) {
        if (group[i] == -1) {
            subgraph.nodes[k].emplace(i);
            subgraph.candidate[i] = 1;
        }
    }

    for (auto i : reverseGraph[s]) {
        if (group[i] == -1) {
            subgraph.nodes[k].emplace(i);
            subgraph.candidate[i] = 1;
        }
    }

    for (auto i : graph[t]) {
        if (group[i] == -1) {
            subgraph.nodes[k].emplace(i);
            subgraph.candidate[i] = 1;
        }
    }

    for (auto i : reverseGraph[t]) {
        if (group[i] == -1) {
            subgraph.nodes[k].emplace(i);
            subgraph.candidate[i] = 1;
        }
    }

    // extract edge
    for (auto i : subgraph.nodes[k]) {
        subgraph.graph[i].clear();
        for (auto v : graph[i]) {
            if (subgraph.candidate[v] == 1) {
                subgraph.graph[i].emplace_back(v);
            }
        }
    }

    for (auto i : subgraph.nodes[k]) {
        subgraph.degree[k][i] = subgraph.graph[i].size();
    }

    for (auto i : subgraph.nodes[k]) {
        subgraph.candidate[i] = 0;
    }
}

void DynamicGraph::FindAllFreeClique(KSubgraph& subgraph, int k, vector<int>& prev, set<Clique*>& freeCliqueList) {
    auto& nodes = subgraph.nodes;
    auto& degree = subgraph.degree;
    auto& graph = subgraph.graph;

    if (k == 2) {
        for (auto i : nodes[k]) {
            for (int j=0;j<degree[k][i];j++) {
                Clique* c = new Clique();
                c->nodes = prev;
                c->nodes.push_back(i);
                c->nodes.push_back(graph[i][j]);
                freeCliqueList.emplace(c);
            }
        }
    } else {
        for (auto i : nodes[k]) {
            if (degree[k][i] < k-1) {
                continue;
            }

            BuildSubGraph(subgraph, k, i);
            prev.push_back(i);
            FindAllFreeClique(subgraph, k-1, prev, freeCliqueList);
            prev.pop_back();
        }
    }
}

void DynamicGraph::MakeSubgraphByFreeEdge(KSubgraph& subgraph, int k, int s, int t) {
    subgraph.nodes[k].clear();

    // mark
    subgraph.nodes[k].emplace(s);
    subgraph.candidate[s] = 1;
    subgraph.nodes[k].emplace(t);
    subgraph.candidate[t] = 1;

    set<int> candS = {s};

    for (auto i : graph[s]) {
        candS.emplace(i);
    }

    for (auto i : reverseGraph[s]) {
        candS.emplace(i);
    }

    for (auto i : graph[t]) {
        if (candS.find(i) != candS.end()) {
            subgraph.nodes[k].emplace(i);
            subgraph.candidate[i] = 1;
        }
    }

    for (auto i : reverseGraph[t]) {
        if (candS.find(i) != candS.end()) {
            subgraph.nodes[k].emplace(i);
            subgraph.candidate[i] = 1;
        }
    }

    // extract edge
    for (auto i : subgraph.nodes[k]) {
        subgraph.graph[i].clear();
        for (auto v : graph[i]) {
            if (subgraph.candidate[v] == 1) {
                subgraph.graph[i].emplace_back(v);
            }
        }
    }

    for (auto i : subgraph.nodes[k]) {
        subgraph.degree[k][i] = subgraph.graph[i].size();
    }

    for (auto i : subgraph.nodes[k]) {
        subgraph.candidate[i] = 0;
    }
}

void DynamicGraph::FindAllCliqueByFreeEdge(KSubgraph& subgraph, int k, vector<int>& prev, set<int>& modifiedGroup, int s1, int t1) {
    auto& nodes = subgraph.nodes;
    auto& degree = subgraph.degree;
    auto& graph = subgraph.graph;

    if (k == 2) {
        for (auto i : nodes[k]) {
            for (int j=0;j<degree[k][i];j++) {
                vector<int> tmp = prev;
                tmp.push_back(i);
                tmp.push_back(graph[i][j]);

                set<int> nodeSet(tmp.begin(), tmp.end());
                if (nodeSet.find(s1)!=nodeSet.end() && nodeSet.find(t1)!=nodeSet.end()) {
                    set<int> type;
                    for (auto p : nodeSet) {
                        type.emplace(group[p]);
                    }
                    if (type.size() == 2) {
                        Clique* c = new Clique();
                        c->nodes = move(tmp);

                        for (auto s : c->nodes) {
                            index[s].emplace(c);
                        }

                        for (auto s : c->nodes) {
                            if (group[s] != -1) {
                                modifiedGroup.emplace(group[s]);
                                break;
                            }
                        }
                    }
                }
            }
        }
    } else {
        for (auto i : nodes[k]) {
            if (degree[k][i] < k-1) {
                continue;
            }

            BuildSubGraph(subgraph, k, i);
            prev.push_back(i);
            FindAllCliqueByFreeEdge(subgraph, k-1, prev, modifiedGroup, s1, t1);
            prev.pop_back();
        }
    }
}

void DynamicGraph::Deletion(int u, int v) {
    // update graph
    DeleteEdge(u, v);

    // rm invalid clique from index
    vector<Clique*> deletedClique;

    for (auto c : index[u]) {
        for (auto node : c->nodes) {
            if (node == v) {
                deletedClique.push_back(c);
                break;
            }
        }
    }

    for (auto c : deletedClique) {
        DeleteIndex(c);
    }

    // if u and v from the same solution, then do something
    if (group[u] == group[v] && group[u] != -1) {
        // rm the invalid clique
        vector<int> nodes = move(invGroup[group[u]]);
        invGroup.erase(group[u]);

        // find candidate cliques
        set<Clique*> cs;
        for (auto i : nodes) {
            group[i] = -1;
            cs.insert(index[i].begin(), index[i].end());
        }

        if (!cs.empty()) {
            // find new cliques
            auto result = DisjointClique(cs);
            vector<int> groupID;

            for (auto targetC : result) {
                for (auto node : targetC) {
                    group[node] = targetC.front();
                }
                invGroup[targetC.front()] = targetC;
                groupID.emplace_back(targetC.front());
            }

            // delete invalid cliques from index
            for (int newGroupID : groupID) {
                for (auto node : invGroup[newGroupID]) {
                    auto nodeIndex = move(index[node]);
                    for (auto i : nodeIndex) {
                        DeleteIndex(i);
                    }
                }
            }

            // update index using new cliques
            for (int newGroupID : groupID) {
                FindIndexByClique(invGroup[newGroupID]);
            }

            // update index by new free nodes
            // collect new free nodes
            set<int> newFreeNodes(nodes.begin(), nodes.end());
            for (auto targetC : result) {
                for (auto node : targetC) {
                    newFreeNodes.erase(node);
                }
            }

            set<int> modifiedGroup(groupID.begin(), groupID.end());

            ExtendIndexByFreeNodes(newFreeNodes, modifiedGroup, modifiedGroup);
            Swap(modifiedGroup);
        } else {
            // update index by new free nodes
            set<int> newFreeNodes(nodes.begin(), nodes.end());
            set<int> modifiedGroup;

            ExtendIndexByFreeNodes(newFreeNodes, modifiedGroup, modifiedGroup);
            Swap(modifiedGroup);
        }
    }
}

void DynamicGraph::DeleteEdge(int u, int v) {
    if (vertexRank[u] < vertexRank[v]) {
        swap(u, v);
    }

    for (int i=0;i<graph[u].size();i++) {
        if (graph[u][i] == v) {
            graph[u].erase(graph[u].begin()+i);
            break;
        }
    }

    for (int i=0;i<reverseGraph[v].size();i++) {
        if (reverseGraph[v][i] == u) {
            reverseGraph[v].erase(reverseGraph[v].begin()+i);
            break;
        }
    } 
}

void DynamicGraph::MakeSubgraphByFreeNodes(KSubgraph& subgraph, int k, set<int>& freeNodes, set<int>& newGroupIDSet) {
    subgraph.nodes[k].clear();

    // mark all one-hop neighbors
    for (auto i : freeNodes) {
        subgraph.nodes[k].emplace(i);
        subgraph.candidate[i] = 1;

        for (auto j : graph[i]) {
            if (newGroupIDSet.find(group[j]) == newGroupIDSet.end()) { // rm redudant computation
                subgraph.nodes[k].emplace(j);
                subgraph.candidate[j] = 1;
            }
        }

        for (auto j : reverseGraph[i]) {
            if (newGroupIDSet.find(group[j]) == newGroupIDSet.end()) {
                subgraph.nodes[k].emplace(j);
                subgraph.candidate[j] = 1;
            }
        }
    }

    // extract edge
    for (auto i : subgraph.nodes[k]) {
        subgraph.graph[i].clear();

        if (group[i] != -1) { // only add group-free or group-same-group edge
            for (auto v : graph[i]) {
                if (subgraph.candidate[v] == 1 && (group[v] == -1 || group[v] == group[i])) {
                    subgraph.graph[i].emplace_back(v);
                }
            }
        } else { // free nodes, add its all edges
            for (auto v : graph[i]) {
                if (subgraph.candidate[v] == 1) {
                    subgraph.graph[i].emplace_back(v);
                }
            }
        }
    }

    for (auto i : subgraph.nodes[k]) {
        subgraph.degree[k][i] = subgraph.graph[i].size();
    }

    for (auto i : subgraph.nodes[k]) {
        subgraph.candidate[i] = 0;
    }
}
void DynamicGraph::FindAllCliqueByFreeNodes(KSubgraph& subgraph, int k, vector<int>& prev, set<int>& modifiedGroup, set<int>& freeNodes) {
    auto& nodes = subgraph.nodes;
    auto& degree = subgraph.degree;
    auto& graph = subgraph.graph;

    if (k == 2) {
        for (auto i : nodes[k]) {
            for (int j=0;j<degree[k][i];j++) {
                vector<int> tmp = prev;
                tmp.push_back(i);
                tmp.push_back(graph[i][j]);

                set<int> nodeSet(tmp.begin(), tmp.end());
                bool containFlag = false;

                // should contain one of the target freeNodes
                {
                    auto it1 = nodeSet.begin();
                    auto it2 = freeNodes.begin();

                    while (it1!=nodeSet.end() && it2!=freeNodes.end()) {
                        if (*it1 == *it2) {
                            containFlag = true;
                            break;
                        } else if (*it1 < *it2) {
                            it1++;
                        } else {
                            it2++;
                        }
                    }
                }

                if (containFlag) {
                    set<int> type;
                    for (auto p : nodeSet) {
                        type.emplace(group[p]);
                    }
                    if (type.size() == 2) {
                        Clique* c = new Clique();
                        c->nodes = move(tmp);

                        for (auto s : c->nodes) {
                            index[s].emplace(c);
                        }

                        for (auto s : c->nodes) {
                            if (group[s] != -1) {
                                modifiedGroup.emplace(group[s]);
                                break;
                            }
                        }
                    }
                }
            }
        }
    } else {
        for (auto i : nodes[k]) {
            if (degree[k][i] < k-1) {
                continue;
            }

            BuildSubGraph(subgraph, k, i);
            prev.push_back(i);
            FindAllCliqueByFreeNodes(subgraph, k-1, prev, modifiedGroup, freeNodes);
            prev.pop_back();
        }
    }
}

void DynamicGraph::ExtendIndexByFreeNodes(set<int>& newFreeNodes, set<int>& completeGroup, set<int>& modifiedGroup) {
    if (!newFreeNodes.empty()) {
        MakeSubgraphByFreeNodes(subgraph, k, newFreeNodes, completeGroup);

        for (auto i : subgraph.nodes[k]) {
            if (subgraph.degree[k][i] < k-1) {
                continue;
            }

            BuildSubGraph(subgraph, k, i);
            vector<int> prev = {i};
            FindAllCliqueByFreeNodes(subgraph, k-1, prev, modifiedGroup, newFreeNodes);
        }
    }
}

void DynamicGraph::DeleteIndex(Clique* c) {
    for (auto i : c->nodes) {
        index[i].erase(c);
    }
    delete c;
}

vector<vector<int>> DynamicGraph::DisjointClique(set<Clique*>& cs) {
    unordered_map<int, int> score;

    for (auto c : cs) {
        for (auto node : c->nodes) {
            score[node]++;
        }
    }

    vector<pair<int, Clique*>> scoreClique;
    for (auto c : cs) {
        int s = 0;
        for (auto i : c->nodes) {
            s += score[i];
        }
        scoreClique.emplace_back(s, c);
    }
    sort(scoreClique.begin(), scoreClique.end());

    for (auto it=score.begin();it!=score.end();it++) {
        it->second = 0;
    }

    vector<vector<int>> result;

    for (auto [s, c] : scoreClique) {
        bool flag = true;

        for (auto i : c->nodes) {
            if (score[i] == 1) {
                flag = false;
                break;
            }
        }

        if (flag) {
            for (auto i : c->nodes) {
                score[i] = 1;
            }
            result.push_back(c->nodes);
        }
    }

    return result;
}

void DynamicGraph::Info() {
    int sum = 0;
    for (auto& i : index) {
        sum += i.size();
    }

    cout << "Index: " << sum / k << endl;
    cout << "Solution: " << invGroup.size() << endl << endl;
}

void DynamicGraph::IndexCorrectness() {
    vector<vector<int>> cliqueList;
    cliqueList.reserve(invGroup.size());

    for (auto [k, c] : invGroup) {
        cliqueList.emplace_back(c);
    }

    for (int i=0;i<n;i++) {
        index[i].clear();
    }

    for (auto c : cliqueList) {
        FindIndexByClique(c);
    }
}