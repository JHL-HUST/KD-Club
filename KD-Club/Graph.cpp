#include "Graph.h"

using namespace std;

void Graph::swap_pos(ui i, ui j)
{
    std::swap(SR[i], SR[j]);
    SR_rid[SR[i]] = i;
    SR_rid[SR[j]] = j;
}

Graph::Graph(const char *_filename, const int _K)
{
    filename = string(_filename);
    fn = filename.substr(filename.rfind('/') + 1);
    K = _K;
    n = m = 0;
    k_defective.clear();

    pstart = nullptr;
    pend = nullptr;
    pend2 = nullptr;
    edges = nullptr;
    edges2 = nullptr;
    exist = nullptr;
    degree = nullptr;
    tri_num = nullptr;
    in_clique = nullptr;

    node_trans = nullptr;
    node_trans_rid = nullptr;
    matrix = nullptr;
    cn = nullptr;
    degree_in_S = nullptr;
    SR = nullptr;
    SR_rid = nullptr;
    neighbors = nullptr;
    levels = nullptr;
    g_color = nullptr;
    node_num = 0;
    branch_counting = 0;
}

Graph::~Graph()
{
    if (pstart != nullptr)
    {
        delete[] pstart;
        pstart = nullptr;
    }
    if (pend != nullptr)
    {
        delete[] pend;
        pend = nullptr;
    }
    if (pend2 != nullptr)
    {
        delete[] pend2;
        pend2 = nullptr;
    }
    if (edges != nullptr)
    {
        delete[] edges;
        edges = nullptr;
    }
    if (edges2 != nullptr)
    {
        delete[] edges2;
        edges2 = nullptr;
    }
    if (exist != nullptr)
    {
        delete[] exist;
        exist = nullptr;
    }
    if (degree != nullptr)
    {
        delete[] degree;
        degree = nullptr;
    }
    if (tri_num != nullptr)
    {
        delete[] tri_num;
        tri_num = nullptr;
    }
    if (matrix != nullptr)
    {
        delete[] matrix;
        matrix = nullptr;
    }
    if (cn != nullptr)
    {
        delete[] cn;
        cn = nullptr;
    }
    if (SR != nullptr)
    {
        delete[] SR;
        SR = nullptr;
    }
    if (SR_rid != nullptr)
    {
        delete[] SR_rid;
        SR_rid = nullptr;
    }
    if (degree_in_S != nullptr)
    {
        delete[] degree_in_S;
        degree_in_S = nullptr;
    }
    if (neighbors != nullptr)
    {
        delete[] neighbors;
        neighbors = nullptr;
    }
    if (node_trans != nullptr)
    {
        delete[] node_trans;
        node_trans = nullptr;
    }
    if (node_trans_rid != nullptr)
    {
        delete[] node_trans_rid;
        node_trans_rid = nullptr;
    }
    if (levels != nullptr)
    {
        delete[] levels;
        levels = nullptr;
    }
    if (g_color != nullptr)
    {
        delete[] g_color;
        g_color = nullptr;
    }
    if (in_clique != nullptr)
    {
        delete[] in_clique;
        in_clique = nullptr;
    }
}

void Graph::read_graph()
{
    FILE *f = Utility::open_file(filename.c_str(), "r");
    char temp;
    for (int i = 0; i < 5; i++)
        fscanf(f, "%1s", &temp);

    fscanf(f, "%u%lu", &n, &m);
    m *= 2;
    printf("*\tn = %s; m = %s (undirected)\n", Utility::integer_to_string(n).c_str(), Utility::integer_to_string(m / 2).c_str());

    vp.clear();
    for (ui i = 0; i < m / 2; i++)
    {
        ui a, b;
        fscanf(f, "%1s%u%u", &temp, &a, &b);
        a--;
        b--;
        if (a >= n || b >= n)
        {
            printf("!!! Vertex IDs must be between 0 and n-1. Exit !!!\n");
            return;
        }
        vp.pb(mp(a, b));
        vp.pb(mp(b, a));
    }
    sort(vp.begin(), vp.end());

    if (pstart != nullptr)
        delete[] pstart;
    pstart = new ept[n + 1];
    if (edges != nullptr)
        delete[] edges;
    if (pend != nullptr)
        delete[] pend;
    pend = new ept[n];
    pend2 = new ept[n];
    edges = new ui[m];
    edges2 = new ui[m];
    tri_num = new ui[m];
    exist = new bool[n];
    in_clique = new bool[n];
    degree = new ui[n];
    g_color = new ui[n];
    rem_node_num = n;
    pstart[0] = 0;
    ui idx = 0;
    for (ui i = 0; i < n; i++)
    {
        pstart[i + 1] = pstart[i];
        while (idx < vp.size() && vp[idx].first == i)
        {
            edges[pstart[i + 1]] = vp[idx].second;
            edges2[pstart[i + 1]] = vp[idx].second;
            idx++;
            pstart[i + 1]++;
        }
        degree[i] = pstart[i + 1] - pstart[i];
    }
    for (ui i = 0; i < n; i++)
    {
        pend[i] = pstart[i + 1];
        pend2[i] = pstart[i + 1];
    }
    fclose(f);

#ifndef NDEBUG
    printf("Finished reading graph\n");
#endif
}

ui Graph::degree_based_heu(ui *nodes, bool *in_cand, ui &s_end, ui &r_end)
{
    if (s_end == r_end)
    {
        return s_end;
    }
    ui max_degree = 0;
    int node_pos = -1;
    for (ui i = s_end; i < r_end; i++)
    {
        ui cur_node = nodes[i];
        int new_degree = 0;
        for (ui j = pstart[cur_node]; j < pstart[cur_node + 1]; j++)
        {
            if (in_cand[edges[j]])
            {
                new_degree++;
            }
        }
        if (new_degree > max_degree || (new_degree == max_degree && node_pos == -1))
        {
            max_degree = new_degree;
            node_pos = i;
        }
    }
    assert(node_pos != -1);
    ui node = nodes[node_pos];
    nodes[node_pos] = nodes[s_end];
    nodes[s_end++] = node;
    r_end = s_end;
    for (ui j = pstart[node]; j < pstart[node + 1]; j++)
    {
        if (in_cand[edges[j]])
        {
            nodes[r_end++] = edges[j];
        }
    }
    memset(in_cand, false, sizeof(bool) * n);
    for (ui j = s_end; j < r_end; j++)
    {
        in_cand[nodes[j]] = true;
    }
    return degree_based_heu(nodes, in_cand, s_end, r_end);
}
void Graph::fastLB1()
{
    cout << "Start FastLB online to calculate a lower bound!" << endl;
    Timer t;
    ui *nodes = new ui[n];
    bool *in_cand = exist;
    ui max_clique = k_defective.size();
    for (ui i = 0; i < n; i++)
    {
        ui degree = pstart[i + 1] - pstart[i];
        if (degree + 1 <= max_clique)
        {
            continue;
        }
        memset(in_cand, false, sizeof(bool) * n);
        nodes[0] = i;
        ui s_end = 1, r_end = 1;
        for (ui j = pstart[i]; j < pstart[i + 1]; j++)
        {
            if (pstart[edges[j] + 1] - pstart[edges[j]] + 1 <= max_clique)
            {
                continue;
            }
            in_cand[edges[j]] = true;
            nodes[r_end++] = edges[j];
        }
        ui cli_size = degree_based_heu(nodes, in_cand, s_end, r_end);
        if (cli_size > max_clique)
        {
            k_defective.clear();
            memset(in_clique, false, sizeof(bool) * n);
            for (ui i = 0; i < s_end; i++)
            {
                k_defective.push_back(nodes[i]);
                in_clique[nodes[i]] = true;
            }
            max_clique = cli_size;
        }
    }
    if (!clique_assert())
    {
        cout << "verify clique failed, please check the codes!" << endl;
    }
    cout << "FastLB online finishes! Find a clique of size " << k_defective.size() << " at time " << Utility::integer_to_string(t.elapsed()).c_str() << "." << endl;
    delete[] nodes;
}
void Graph::fastLB2()
{
    cout << "Start FastLB offline to calculate a lower bound!" << endl;
    Timer t;
    ui *nodes = new ui[n];
    bool *in_cand = exist;
    ui max_clique = 0;
    for (ui i = 0; i < n; i++)
    {
        ui degree = pstart[i + 1] - pstart[i];
        if (degree + 1 <= max_clique)
        {
            continue;
        }
        memset(in_cand, false, sizeof(bool) * n);
        nodes[0] = i;
        ui s_end = 1, r_end = 1;
        for (ui j = pstart[i]; j < pstart[i + 1]; j++)
        {
            if (pstart[edges[j] + 1] - pstart[edges[j]] + 1 <= max_clique)
            {
                continue;
            }
            in_cand[edges[j]] = true;
            nodes[r_end++] = edges[j];
        }
        ui cli_size = degree_based_heu2(nodes, in_cand, s_end, r_end);
        if (cli_size > max_clique)
        {
            k_defective.clear();
            memset(in_clique, false, sizeof(bool) * n);
            for (ui i = 0; i < s_end; i++)
            {
                k_defective.push_back(nodes[i]);
                in_clique[nodes[i]] = true;
            }
            max_clique = cli_size;
        }
    }
    if (!clique_assert())
    {
        cout << "verify clique failed, please check the codes!" << endl;
    }
    cout << "FastLB offline finishes! Find a clique of size " << k_defective.size() << " at time " << Utility::integer_to_string(t.elapsed()).c_str() << "." << endl;
    delete[] nodes;
}
ui Graph::degree_based_heu2(ui *nodes, bool *in_cand, ui &s_end, ui &r_end)
{
    if (s_end == r_end)
    {
        return s_end;
    }
    ui max_degree = 0;
    int node_pos = -1;
    for (ui i = s_end; i < r_end; i++)
    {
        ui cur_node = nodes[i];
        ui degree = pstart[cur_node + 1] - pstart[cur_node];
        if (degree > max_degree || node_pos == -1)
        {
            max_degree = degree;
            node_pos = i;
        }
    }
    assert(node_pos != -1);
    ui node = nodes[node_pos];
    nodes[node_pos] = nodes[s_end];
    nodes[s_end] = node;
    s_end++;
    r_end = s_end;
    for (ui j = pstart[node]; j < pstart[node + 1]; j++)
    {
        if (in_cand[edges[j]])
        {
            nodes[r_end++] = edges[j];
        }
    }
    memset(in_cand, false, sizeof(bool) * n);
    for (ui j = s_end; j < r_end; j++)
    {
        in_cand[nodes[j]] = true;
    }
    return degree_based_heu2(nodes, in_cand, s_end, r_end);
}

bool Graph::clique_assert()
{
    bool *in_cand = exist;
    memset(in_cand, false, sizeof(bool) * n);
    for (ui i = 0; i < k_defective.size(); i++)
    {
        for (ui j = pstart[k_defective[i]]; j < pstart[k_defective[i] + 1]; j++)
        {
            exist[edges[j]] = true;
        }
        for (ui j = i + 1; j < k_defective.size(); j++)
        {
            if (!exist[k_defective[j]])
            {
                return false;
            }
        }
        for (ui j = pstart[k_defective[i]]; j < pstart[k_defective[i] + 1]; j++)
        {
            exist[edges[j]] = false;
        }
    }
    return true;
}

bool Graph::can_color(ui u, std::set<ui> &c, bool *in_neigh, int mode)
{
    for (ui i = pstart[u]; i < pend2[u]; i++)
    {
        ui v = edges2[i];
        if ((!in_neigh[v] && mode == 1) || (mode == 2 && in_neigh[v]))
        {
            continue;
        }
        if (c.find(v) != c.end())
        {
            return false;
        }
    }
    return true;
}
ui Graph::cal_col_bound_node(ui u)
{
    bool *is_neighbor = new bool[n];
    memset(is_neighbor, false, sizeof(bool) * n);
    vector<set<ui>> colors;
    int color_num = 0;
    colors.push_back(set<ui>());
    for (ui i = pstart[u]; i < pend2[u]; i++)
    {
        ui v = edges2[i];
        if (!exist[v])
        {
            continue;
        }
        is_neighbor[v] = true;
    }
    for (ui i = pstart[u]; i < pend2[u]; i++)
    {
        ui v = edges2[i];
        if (!exist[v])
        {
            continue;
        }
        ui cur_color = 0;
        for (; cur_color < color_num; cur_color++)
        {
            if (can_color(v, colors[cur_color], is_neighbor, 1))
            {
                colors[cur_color].insert(v);
                break;
            }
        }
        if (cur_color == color_num)
        {
            colors.push_back(set<ui>());
            color_num++;
            colors[cur_color].insert(v);
        }
    }
    ui cur_ub = 1 + color_num;
    map<ui, ui> addk_num_pair;
    int remain = degree[u] - color_num;
    if (min(color_num, remain) >= K)
    {
        delete[] is_neighbor;
        return cur_ub + K;
    }
    if (remain >= K && cur_ub + K <= k_defective.size())
    {
        delete[] is_neighbor;
        return cur_ub + K;
    }
    ui addk = 1;
    while (remain > 0)
    {
        if (addk_num_pair.find(addk) == addk_num_pair.end())
        {
            addk_num_pair[addk] = min(color_num, remain);
        }
        else
        {
            addk_num_pair[addk] += min(color_num, remain);
        }
        remain -= color_num;
        addk++;
    }
    is_neighbor[u] = true;

    colors.clear();
    colors.push_back(set<ui>());
    color_num = 0;
    ui counter = 0;
    for (ui i = 0; i < n; i++)
    {
        if (is_neighbor[i] || !exist[i])
        {
            continue;
        }
        ui c = 0;
        counter++;
        for (; c < color_num; c++)
        {
            if (can_color(i, colors[c], is_neighbor, 2))
            {
                colors[c].insert(i);
                break;
            }
        }
        if (c == color_num)
        {
            colors.push_back(set<ui>());
            color_num++;
            colors[c].insert(i);
        }
    }
    addk = 1;
    remain = counter;
    while (remain > 0)
    {
        if (addk_num_pair.find(addk) == addk_num_pair.end())
        {
            addk_num_pair[addk] = min(color_num, remain);
        }
        else
        {
            addk_num_pair[addk] += min(color_num, remain);
        }
        remain -= color_num;
        addk++;
    }
    ui cur_k = K;
    for (ui i = 1; i <= K; i++)
    {
        if (addk_num_pair.find(i) == addk_num_pair.end())
            continue;
        if (addk_num_pair[i] * i >= cur_k)
        {
            cur_ub += cur_k / i;
            break;
        }
        else
        {
            cur_k -= addk_num_pair[i] * i;
            cur_ub += addk_num_pair[i];
        }
    }
    delete[] is_neighbor;
    return cur_ub;
}
void Graph::check_vertex(ui *Qv, ui &Qv_n, ui lb, bool &removed, bool *marked_for_node, bool *marked_for_edge, bool first_time)
{
    cout << "  Start checking useless nodes!" << endl;
    removed = false;
    ui removed_node_count = 0;
    bool *in_queue = new bool[n];
    memset(in_queue, false, sizeof(bool) * n);
    bool *need_check = new bool[n];
    memset(need_check, false, sizeof(bool) * n);
    for (ui u = 0; u < n; u++)
    {
        if (!exist[u] || !marked_for_node[u] || in_clique[u])
        {
            continue;
        }
        ui ub1 = 1 + degree[u] + min(rem_node_num - degree[u] - 1, K);
        if (ub1 > lb && !first_time)
        {
            need_check[u] = true;
        }
        if (ub1 <= lb)
        {
            exist[u] = false;
            removed_node_count++;
            rem_node_num--;
            for (ui j = pstart[u]; j < pend2[u]; j++)
            {
                ui v = edges2[j];
                if (!exist[v] || in_clique[v])
                {
                    continue;
                }
                marked_for_edge[v] = true;
                degree[v]--;
                ui ub_v = 1 + degree[v] + min(rem_node_num - degree[v] - 1, K);
                if (v < u && !in_queue[v] && ub_v <= lb)
                {
                    Qv[Qv_n++] = v;
                    in_queue[v] = true;
                }
                else if (v < u && !in_queue[v] && ub_v > lb)
                {
                    need_check[v] = true;
                }
            }
        }
    }
    for (ui u = 0; u < n; u++)
    {
        if (!exist[u] || !need_check[u] || in_queue[u] || in_clique[u])
        {
            continue;
        }
        ui ub1 = 1 + degree[u] + min(rem_node_num - degree[u] - 1, K);
        ui ub2 = ub1;
        if (ub1 >= lb && !first_time)
            ub2 = cal_col_bound_node(u);
        if (ub1 < lb || ub2 < lb)
        {
            Qv[Qv_n++] = u;
            in_queue[u] = true;
        }
    }
    memset(need_check, false, sizeof(bool) * n);
    for (ui i = 0; i < Qv_n; i++)
    {
        ui u = Qv[i];
        assert(in_queue[u] && exist[u] && !in_clique[u]);
        exist[u] = false;
        removed_node_count++;
        rem_node_num--;
        for (ui j = pstart[u]; j < pend2[u]; j++)
        {
            ui v = edges2[j];
            if (!exist[v] || in_clique[v])
            {
                continue;
            }
            marked_for_edge[v] = true;
            degree[v]--;
            ui ub_v1 = 1 + degree[v] + min(rem_node_num - degree[v] - 1, K);
            if (!in_queue[v] && ub_v1 <= lb)
            {
                Qv[Qv_n++] = v;
                in_queue[v] = true;
            }
            else if (!in_queue[v] && ub_v1 > lb)
            {
                need_check[v] = true;
            }
        }
        in_queue[u] = false;
        if (first_time)
            continue;
        if (i == Qv_n - 1)
        {
            for (ui i = 0; i < n; i++)
            {
                if (!need_check[i] || !exist[i])
                {
                    continue;
                }
                ui ub1 = 1 + degree[i] + min(rem_node_num - degree[i] - 1, K);
                ui ub2 = ub1;
                if (ub1 >= lb)
                    ub2 = cal_col_bound_node(i);
                if (ub1 <= lb || ub2 <= lb)
                {
                    Qv[Qv_n++] = i;
                    in_queue[i] = true;
                }
            }
            memset(need_check, false, sizeof(bool) * n);
        }
    }
    delete[] need_check;
    for (ui u = 0; u < n; u++)
    {
        if (!exist[u])
        {
            degree[u] = 0;
            pend2[u] = pstart[u];
            continue;
        }
        ui temp_end = pend2[u];
        pend2[u] = pstart[u];
        for (ui j = pstart[u]; j < temp_end; j++)
        {
            ui v = edges2[j];
            if (!exist[v])
            {
                continue;
            }
            if (!first_time)
            {
                if (u < v)
                {
                    edges2[pend2[u]++] = v;
                }
            }
            else
            {
                edges2[pend2[u]++] = v;
            }
        }
    }
    delete[] in_queue;
    if (removed_node_count > 0)
    {
        removed = true;
    }
    cout << "    After checking vertex, " << removed_node_count << " nodes are removed! " << rem_node_num << " nodes remain!" << endl;
}
ui Graph::cal_col_bound_edge(ui u, ui v, ui edge_pos)
{
    bool *both_neigh = new bool[n];
    memset(both_neigh, false, sizeof(bool) * n);
    vector<ui> both_neigh_nodes;
    bool *one_neigh = new bool[n];
    memset(one_neigh, false, sizeof(bool) * n);
    vector<ui> one_neigh_nodes;
    for (ui i = pstart[u]; i < pend2[u]; i++)
    {
        if (!exist[edges2[i]] || edges2[i] == v)
        {
            continue;
        }
        one_neigh[edges2[i]] = true;
    }
    for (ui i = pstart[v]; i < pend2[v]; i++)
    {
        if (!exist[edges2[i]] || edges2[i] == u)
        {
            continue;
        }
        if (one_neigh[edges2[i]])
        {
            both_neigh[edges2[i]] = true;
            one_neigh[edges2[i]] = false;
            both_neigh_nodes.push_back(edges2[i]);
        }
        else
        {
            one_neigh[edges2[i]] = true;
            one_neigh_nodes.push_back(edges2[i]);
        }
    }
    for (ui i = pstart[u]; i < pend2[u]; i++)
    {
        if (!exist[edges2[i]] || edges2[i] == v)
        {
            continue;
        }
        if (one_neigh[edges2[i]])
        {
            one_neigh_nodes.push_back(edges2[i]);
        }
        else
        {
            assert(both_neigh[edges2[i]]);
        }
    }
    map<ui, ui> addk_num_pair;
    vector<set<ui>> colors;
    colors.clear();
    colors.push_back(set<ui>());
    int color_num = 0;
    for (ui node : both_neigh_nodes)
    {
        ui cur_color = 0;
        for (; cur_color < color_num; cur_color++)
        {
            if (can_color(node, colors[cur_color], both_neigh, 1))
            {
                colors[cur_color].insert(node);
                break;
            }
        }
        if (cur_color == color_num)
        {
            color_num++;
            colors.push_back(set<ui>());
            colors[cur_color].insert(node);
        }
    }

    assert(both_neigh_nodes.size() == tri_num[edge_pos]);
    if (both_neigh_nodes.size() - color_num >= K)
    {
        delete[] both_neigh;
        delete[] one_neigh;
        return 2 + color_num + K;
    }
    ui add_k = 0;
    int remain = both_neigh_nodes.size();
    while (remain > 0)
    {
        if (addk_num_pair.find(add_k) == addk_num_pair.end())
        {
            addk_num_pair[add_k] = min(remain, color_num);
        }
        else
        {
            addk_num_pair[add_k] += min(remain, color_num);
        }
        remain -= color_num;
        add_k++;
    }
    ui cur_ub = 2 + both_neigh_nodes.size();
    ui cur_k = K - (both_neigh_nodes.size() - color_num);
    colors.clear();
    colors.push_back(set<ui>());
    color_num = 0;
    for (ui node : one_neigh_nodes)
    {
        ui cur_color = 0;
        for (; cur_color < color_num; cur_color++)
        {
            if (can_color(node, colors[cur_color], one_neigh, 1))
            {
                colors[cur_color].insert(node);
                break;
            }
        }
        if (cur_color == color_num)
        {
            color_num++;
            colors.push_back(set<ui>());
            colors[cur_color].insert(node);
        }
    }
    if (color_num >= cur_k)
    {
        delete[] both_neigh;
        delete[] one_neigh;
        return cur_ub + cur_k;
    }
    cur_k -= color_num;
    cur_ub += color_num;
    if (2 * (one_neigh_nodes.size() - color_num) >= cur_k)
    {
        delete[] both_neigh;
        delete[] one_neigh;
        return cur_ub + cur_k / 2;
    }
    add_k = 1;
    remain = one_neigh_nodes.size();
    while (remain > 0)
    {
        if (addk_num_pair.find(add_k) == addk_num_pair.end())
        {
            addk_num_pair[add_k] = min(remain, color_num);
        }
        else
        {
            addk_num_pair[add_k] += min(remain, color_num);
        }
        remain -= color_num;
        add_k++;
    }
    cur_ub += one_neigh_nodes.size() - color_num;
    cur_k -= (2 * (one_neigh_nodes.size() - color_num));
    colors.clear();
    colors.push_back(set<ui>());
    color_num = 0;
    ui counter = 0;
    one_neigh[u] = true;
    one_neigh[v] = true;
    for (ui i = 0; i < n; i++)
    {
        if (both_neigh[i])
        {
            one_neigh[i] = true;
        }
        if (exist[i] && !one_neigh[i])
        {
            counter++;
        }
    }
    ui debug_count = 0;
    for (ui i = 0; i < n; i++)
    {
        if (!exist[i] || one_neigh[i])
        {
            continue;
        }
        debug_count++;
        ui cur_color = 0;
        for (; cur_color < color_num; cur_color++)
        {
            if (can_color(i, colors[cur_color], one_neigh, 2))
            {
                colors[cur_color].insert(i);
                break;
            }
        }
        if (cur_color == color_num)
        {
            color_num++;
            colors.push_back(set<ui>());
            colors[cur_color].insert(i);
        }
    }
    assert(counter == debug_count);
    delete[] both_neigh;
    delete[] one_neigh;
    add_k = 2;
    remain = counter;
    while (remain > 0)
    {
        if (addk_num_pair.find(add_k) == addk_num_pair.end())
        {
            addk_num_pair[add_k] = min(remain, color_num);
        }
        else
        {
            addk_num_pair[add_k] += min(remain, color_num);
        }
        remain -= color_num;
        add_k++;
    }
    cur_ub = 2;
    if (addk_num_pair.find(0) != addk_num_pair.end())
    {
        cur_ub += addk_num_pair[0];
    }
    cur_k = K;
    for (ui i = 1; i <= K; i++)
    {
        if (addk_num_pair.find(i) == addk_num_pair.end())
            continue;
        if (addk_num_pair[i] * i >= cur_k)
        {
            return cur_ub + cur_k / i;
        }
        else
        {
            cur_k -= addk_num_pair[i] * i;
            cur_ub += addk_num_pair[i];
        }
    }
    return cur_ub;
}
void Graph::triangle_count()
{
    ui *adj = degree;
    memset(tri_num, 0, sizeof(ui) * m);
    memset(adj, 0, sizeof(ui) * n);
    for (ui u = 0; u < n; u++)
    {
        if (!exist[u])
            continue;
        for (ui i = pstart[u]; i < pend2[u]; i++)
        {
            assert(exist[edges2[i]]);
            adj[edges2[i]] = i + 1;
        }
        for (ui i = pstart[u]; i < pend2[u]; i++)
        {
            ui v = edges2[i];
            for (ui j = pstart[v]; j < pend2[v]; j++)
            {
                ui w = edges2[j];
                if (adj[w])
                {
                    tri_num[j]++;
                    tri_num[i]++;
                    tri_num[adj[w] - 1]++;
                }
            }
        }
        for (ui i = pstart[u]; i < pend2[u]; i++)
        {
            adj[edges2[i]] = 0;
        }
    }
    for (ui u = 0; u < n; u++)
    {
        if (!exist[u])
            continue;
        for (ui i = pstart[u]; i < pend2[u]; i++)
        {
            ui v = edges2[i];
            assert(exist[v]);
            if (v > u)
            {
                edges2[pend2[v]] = u;
                tri_num[pend2[v]] = tri_num[i];
                pend2[v]++;
            }
        }
    }
}
void Graph::update_tri_num(ui u, ui v)
{
    ui *adj = degree;
    memset(adj, 0, sizeof(ui) * n);
    for (ui i = pstart[u]; i < pend2[u]; i++)
    {
        adj[edges2[i]] = i + 1;
    }
    for (ui j = pstart[v]; j < pend2[v]; j++)
    {
        if (edges2[j] == u)
        {
            --pend2[v];
            edges2[j] = edges2[pend2[v]];
            tri_num[j] = tri_num[pend2[v]];
        }
        ui w = edges2[j];
        if (adj[w])
        {
            tri_num[j]--;
            tri_num[adj[w] - 1]--;
            for (ui i = pstart[w]; i < pend2[w]; i++)
            {
                if (edges2[i] == u || edges2[i] == v)
                {
                    tri_num[i]--;
                }
            }
        }
    }
}
void Graph::check_edges(ui *Qe, ui &Qe_n, ui lb, bool &removed, bool *marked_for_node, bool *marked_for_edge, bool first_time)
{
    cout << "  Start checking useless edges!" << endl;
    removed = false;
    ui reomved_edge_count = 0;
    bool *in_queue = new bool[n];
    memset(in_queue, false, sizeof(bool) * n);
    for (ui u = 0; u < n; u++)
    {
        if (!exist[u] || !marked_for_edge[u])
            continue;
        for (ui j = pstart[u]; j < pend2[u]; j++)
        {
            ui v = edges2[j];
            if (!exist[v] || v < u || (in_clique[u] && in_clique[v]))
                continue;
            ui ub1 = 2 + tri_num[j] + min(K, rem_node_num - tri_num[j] - 2);
            ui ub2 = ub1;
            if (ub1 > lb && !first_time)
            {
                ub2 = cal_col_bound_edge(u, v, j);
            }
            if (ub1 <= lb || ub2 <= lb)
            {
                removed = true;
                --pend2[u];
                edges2[j] = edges2[pend2[u]];
                tri_num[j] = tri_num[pend2[u]];
                --j;
                reomved_edge_count++;
                update_tri_num(u, v);
                if (!in_queue[u])
                {
                    Qe[Qe_n++] = u;
                    in_queue[u] = true;
                }
                marked_for_node[u] = true;
                if (!in_queue[v])
                {
                    Qe[Qe_n++] = v;
                    in_queue[v] = true;
                }
                marked_for_node[v] = true;
            }
        }
    }
    for (ui i = 0; i < Qe_n; i++)
    {
        ui u = Qe[i];
        assert(exist[u] && in_queue[u]);
        bool remove = false;
        in_queue[u] = false;
        for (ui j = pstart[u]; j < pend2[u]; j++)
        {
            ui v = edges2[j];
            if (!exist[v] || (in_clique[u] && in_clique[v]))
                continue;

            ui ub1 = 2 + tri_num[j] + min(K, rem_node_num - tri_num[j] - 2);
            ui ub2 = ub1;
            if (ub1 > lb && !first_time)
            {
                ub2 = cal_col_bound_edge(u, v, j);
            }
            if (ub1 <= lb || ub2 <= lb)
            {
                --pend2[u];
                edges2[j] = edges2[pend2[u]];
                tri_num[j] = tri_num[pend2[u]];
                --j;
                remove = true;
                reomved_edge_count++;
                update_tri_num(u, v);
                if (!in_queue[v])
                {
                    Qe[Qe_n++] = v;
                    in_queue[v] = true;
                }
                marked_for_node[v] = true;
            }
        }
        if (remove && !in_queue[u])
        {
            Qe[Qe_n++] = u;
            in_queue[u] = true;
        }
    }
    cout << "    After checking edge, " << reomved_edge_count << " edges are removed!" << endl;
    delete[] in_queue;
}
ui Graph::preprocess()
{
    cout << "Preprocessing start!" << endl;
    Timer t;
    bool removed_vertex = false;
    bool removed_edge = false;
    bool first_time = true;
    memset(exist, true, sizeof(bool) * n);
    ui lb = k_defective.size();
    bool *marked_for_node = new bool[n];
    bool *marked_for_edge = new bool[n];
    memset(marked_for_node, true, sizeof(bool) * n);
    memset(marked_for_edge, true, sizeof(bool) * n);
    ui Qv_n = 0;
    ui *Qv = new ui[n];
    check_vertex(Qv, Qv_n, lb, removed_vertex, marked_for_node, marked_for_edge, true);
    delete[] Qv;
    while (first_time || removed_vertex || removed_edge)
    {
        if (!first_time && !removed_edge)
        {
            break;
        }
        ui Qv_n = 0;
        ui *Qv = new ui[n];
        check_vertex(Qv, Qv_n, lb, removed_vertex, marked_for_node, marked_for_edge, false);
        delete[] Qv;
        triangle_count();
        if (!first_time && !removed_vertex)
        {
            break;
        }
        ui Qe_n = 0;
        ui *Qe = new ui[m];
        memset(marked_for_node, false, sizeof(bool) * n);
        check_edges(Qe, Qe_n, lb, removed_edge, marked_for_node, marked_for_edge, first_time);
        delete[] Qe;
        memset(marked_for_edge, false, sizeof(bool) * n);
        memset(degree, 0, sizeof(ui) * n);
        ui counter = 0;
        for (ui u = 0; u < n; u++)
        {
            if (exist[u])
            {
                degree[u] = pend2[u] - pstart[u];
                counter += degree[u];
            }
        }
        first_time = false;
    }
    delete[] marked_for_node;
    delete[] marked_for_edge;
    delete[] in_clique;
    in_clique = nullptr;
    rem_edge_num = 0;
    rem_node_num = 0;
    for (ui u = 0; u < n; u++)
    {
        if (!exist[u])
            continue;
        rem_node_num++;
        for (ui i = pstart[u]; i < pend2[u]; i++)
        {
            if (exist[edges2[i]])
                rem_edge_num++;
        }
    }
    rem_edge_num /= 2;
    pre_time = t.elapsed() * 0.000001;
    cout << "Preprocessing finished! Remaining node num:" << rem_node_num << "; Remaining edge num:" << rem_edge_num << '.' << endl;
    cout << "Preprocessing time:" << pre_time << "s" << endl;
    pre_node = rem_node_num;
    pre_edge = rem_edge_num;
    return construct_bnb_solver(rem_node_num);
}
void Graph::output_pre_result()
{
    string wf = "./file/Realworld/" + fn.substr(0, fn.find('.')) + "-" + to_string(K) + ".col";
    FILE *fout = Utility::open_file(wf.c_str(), "w");
    fprintf(fout, "p edge %u %u\n", node_num, rem_edge_num);
    ui count = 0;
    for (ui i = 0; i < node_num; i++)
    {
        for (ui j = i + 1; j < node_num; j++)
        {
            if (matrix[SR[i] * node_num + SR[j]])
            {
                fprintf(fout, "e %u %u\n", SR[i] + 1, SR[j] + 1);
                count++;
            }
        }
    }
    fclose(fout);
    assert(count == rem_edge_num);
}
ui Graph::construct_bnb_solver(ui rem_node_num)
{
    cout << "Constructing bnb solver!" << endl;
    node_trans = new ui[rem_node_num];
    node_trans_rid = new ui[n];
    SR = new ui[rem_node_num];
    SR_rid = new ui[rem_node_num];
    matrix = new char[rem_node_num * rem_node_num];
    cn = new ui[rem_node_num * rem_node_num];
    degree_in_S = new ui[rem_node_num];
    levels = new ui[rem_node_num];
    neighbors = new ui[rem_node_num];
    g_color = new ui[rem_node_num];
    is_nwoe = new bool[rem_node_num];
    memset(degree_in_S, 0, sizeof(ui) * rem_node_num);
    memset(matrix, 0, sizeof(char) * rem_node_num * rem_node_num);
    memset(cn, 0, sizeof(ui) * rem_node_num * rem_node_num);
    memset(is_nwoe, false, sizeof(bool) * rem_node_num);
    ui R_end = 0;

    cur_uncon_num = 0;
    removed_edge_n = 0;
    for (ui i = 0; i < n; i++)
    {
        if (exist[i])
        {
            // node_trans:ori->node_in_solver
            node_trans[R_end] = i;     //(new->old)
            node_trans_rid[i] = R_end; //(old->new)
            // SR:node_in_solver(use in matrix and cn finding)
            SR[R_end] = R_end;     // pos in matrix
            SR_rid[R_end] = R_end; // name of pos in matrix
            if (degree[i] == 0)
            {
                is_nwoe[R_end] = true;
            }
            R_end++;
        }
    }
    memset(degree, 0, sizeof(ui) * rem_node_num);
    node_num = R_end;
    assert(node_num == rem_node_num);
    for (ui i = 0; i < R_end; i++)
    {
        ui u = node_trans[SR[i]];
        char *matrix_t = matrix + SR[i] * node_num;
        assert(exist[u]);
        for (ui j = pstart[u]; j < pend2[u]; j++)
        {
            ui v = edges2[j];
            if (!exist[v] || v < u)
                continue;
            ui v_ = node_trans_rid[v];
            matrix_t[v_] = 1;
            degree[v_]++;
            matrix[v_ * node_num + SR[i]] = 1;
            degree[SR[i]]++;
        }
    }
    for (ui i = 0; i < R_end; i++)
    {
        ui neighbors_n = 0;
        char *t_matrix = matrix + SR[i] * node_num;
        for (ui j = 0; j < R_end; j++)
        {
            if (t_matrix[SR[j]])
                neighbors[neighbors_n++] = SR[j];
        }
        for (ui j = 0; j < neighbors_n; j++)
            for (ui k = j + 1; k < neighbors_n; k++)
            {
                ++cn[neighbors[j] * node_num + neighbors[k]];
                ++cn[neighbors[k] * node_num + neighbors[j]];
            }
    }
    return R_end;
}
void Graph::show(ui S_end, ui R_end)
{
    cout << "fix part:" << endl;
    for (ui i = 0; i < S_end; i++)
    {
        cout << SR[i] << ' ';
    }
    cout << endl;
    cout << "cand part:" << endl;
    for (ui i = S_end; i < R_end; i++)
    {
        cout << SR[i] << ' ';
    }
    cout << endl;
}

ui Graph::cal_col_bound_reduction(ui u, ui S_end, ui R_end)
{
    char *matrix_t = matrix + u * node_num;
    set<ui> cur_color;
    ui color_num = 0;
    ui counter = 0;
    for (ui i = S_end; i < R_end; i++)
    {
        if (!matrix_t[SR[i]])
        {
            continue;
        }
        counter++;
        if (cur_color.find(g_color[SR[i]]) == cur_color.end())
        {
            cur_color.insert(g_color[SR[i]]);
            color_num++;
        }
    }
    if (counter - color_num >= K - cur_uncon_num)
    {
        return 1 + color_num + K - cur_uncon_num;
    }
    ui cur_ub = 1 + counter;
    ui cur_k = K - cur_uncon_num - (counter - color_num);
    counter = 0;
    color_num = 0;
    cur_color.clear();
    for (ui i = S_end; i < R_end; i++)
    {
        if (matrix_t[SR[i]] || SR[i] == u)
        {
            continue;
        }
        counter++;
        if (cur_color.find(g_color[SR[i]]) == cur_color.end())
        {
            cur_color.insert(g_color[SR[i]]);
            color_num++;
        }
    }
    if (color_num >= cur_k)
    {
        return cur_ub + cur_k;
    }
    cur_k -= color_num;
    cur_ub += color_num;
    return cur_ub + min(counter - color_num, cur_k / 2);
}
void Graph::reduction(ui S_end, ui &R_end, ui last_node, ui level)
{
    int lb = k_defective.size() - S_end;
    ui *Q = new ui[node_num];
    bool *in_queue = new bool[node_num];
    memset(in_queue, 0, sizeof(bool) * node_num);
    ui Q_n = 0;
    char *matrix_t = matrix + last_node * node_num;
    for (ui i = S_end; i < R_end; i++)
    {
        if (matrix_t[SR[i]])
        {
            Q[Q_n++] = SR[i];
            in_queue[SR[i]] = true;
        }
    }
    if (SR_rid[last_node] < S_end)
    {
        ui u = last_node;
        ui *cn_t = cn + u * node_num;
        for (ui j = S_end; j < R_end; j++)
        {
            ui v = SR[j];
            int ub = cn_t[v] + 1 + min(K - cur_uncon_num, R_end - S_end - cn_t[v] - 1);
            if (ub <= lb || S_end - degree_in_S[v] > K - cur_uncon_num)
            {
                swap_pos(j, R_end - 1);
                levels[v] = level;
                R_end--;
                j = j - 1;
                matrix_t = matrix + v * node_num;
                ui neighbors_n = 0;
                for (ui i = S_end; i < R_end; i++)
                {
                    if (matrix_t[SR[i]])
                    {
                        degree[SR[i]]--;
                        neighbors[neighbors_n++] = SR[i];
                        if (!in_queue[SR[i]])
                        {
                            Q[Q_n++] = SR[i];
                            in_queue[SR[i]] = true;
                        }
                    }
                }
                for (ui i = 0; i < neighbors_n; i++)
                {
                    for (ui k = i + 1; k < neighbors_n; k++)
                    {
                        --cn[neighbors[i] * node_num + neighbors[k]];
                        --cn[neighbors[k] * node_num + neighbors[i]];
                    }
                }
            }
        }
    }
    // show(S_end, R_end);
    for (ui i = 0; i < Q_n; i++)
    {
        ui w = Q[i];
        if (SR_rid[w] >= R_end)
        {
            continue;
        }
        ui *cn_t = cn + w * node_num;
        matrix_t = matrix + w * node_num;
        for (ui j = S_end; j < R_end; j++)
        {
            ui v = SR[j];
            if (matrix_t[v] == 0)
                continue;
            int ub = 2 + cn_t[v] + min(K - cur_uncon_num, R_end - S_end - cn_t[v] - 2);
            if (ub <= lb)
            {
                if (removed_edges.size() == removed_edge_n)
                {
                    removed_edges.push_back(std::make_pair(v, w));
                    ++removed_edge_n;
                }
                else
                    removed_edges[removed_edge_n++] = std::make_pair(v, w);
                matrix[v * node_num + w] = matrix[w * node_num + v] = 0;
                degree[v]--;
                degree[w]--;
                ui neighbor_n = 0;

                for (ui k = S_end; k < R_end; k++)
                {
                    if (matrix_t[SR[k]])
                    {
                        ui u = SR[k];
                        cn[u * node_num + v]--;
                        cn[v * node_num + u]--;
                    }
                }
                for (ui k = S_end; k < R_end; k++)
                {
                    if (matrix[node_num * v + SR[k]])
                    {
                        ui u = SR[k];
                        cn[u * node_num + w]--;
                        cn[w * node_num + u]--;
                    }
                }
                ui degree_in_r = degree[w] - degree_in_S[w];
                int ub = 1 + degree_in_r + min(K - cur_uncon_num, R_end - S_end - 1 - degree_in_r);
                if (ub <= lb)
                {
                    swap_pos(SR_rid[w], R_end - 1);
                    R_end--;
                    levels[w] = level;
                    ui neighbors_n = 0;
                    for (ui k = S_end; k < R_end; k++)
                    {
                        if (matrix_t[SR[k]])
                        {
                            neighbors[neighbors_n++] = SR[k];
                            degree[SR[k]]--;
                        }
                    }
                    for (ui p = 0; p < neighbors_n; p++)
                    {
                        for (ui k = p + 1; k < neighbors_n; k++)
                        {
                            --cn[neighbors[p] * node_num + neighbors[k]];
                            --cn[neighbors[k] * node_num + neighbors[p]];
                        }
                    }
                    break;
                }
            }
        }
    }

    delete[] Q;
    delete[] in_queue;
}

ui Graph::cal_color_num(vector<ui> &nodes)
{
    vector<vector<ui>> colors;
    ui color_num = 1;
    colors.push_back(vector<ui>());
    if (nodes.empty())
    {
        return 0;
    }
    for (ui cur_node : nodes)
    {
        ui cur_color = 0;
        for (; cur_color < color_num; cur_color++)
        {
            bool can_color = true;
            for (ui node : colors[cur_color])
            {
                if (matrix[cur_node * node_num + node])
                {
                    can_color = false;
                    break;
                }
            }
            if (can_color)
            {
                colors[cur_color].push_back(cur_node);
                break;
            }
        }
        if (cur_color == color_num)
        {
            color_num++;
            colors.push_back(vector<ui>());
            colors[cur_color].push_back(cur_node);
        }
    }
    return color_num;
}
ui Graph::cal_candidate_bound_colorset(ui S_end, ui R_end)
{
    vector<vector<ui>> partition_by_ucn;
    for (int i = 0; i <= K; i++)
    {
        partition_by_ucn.push_back(vector<ui>());
    }
    ui min_degree = node_num;
    vector<ui> min_degree_node;
    ui count_nwoe = 0;
    for (ui i = S_end; i < R_end; i++)
    {
        if (degree[SR[i]] == 0)
        {
            count_nwoe++;
            continue;
        }
        int uc_num = S_end - degree_in_S[SR[i]];
        assert(uc_num <= K && uc_num >= 0);
        partition_by_ucn[uc_num].push_back(SR[i]);
        ui degree_in_r = degree[SR[i]] - degree_in_S[SR[i]];
        assert(degree_in_r >= 0 && degree_in_r < R_end - S_end);
        if (degree_in_r < min_degree)
        {
            min_degree_node.clear();
            min_degree_node.push_back(i);
            min_degree = degree_in_r;
        }
        else if (degree_in_r == min_degree)
        {
            min_degree_node.push_back(i);
        }
    }
    ui count_num = S_end;
    ui cur_add_uc = 0;
    bool flag = false;
    for (int i = 0; i <= K; i++)
    {
        if (cur_add_uc + partition_by_ucn[i].size() * i <= K - cur_uncon_num)
        {
            cur_add_uc += partition_by_ucn[i].size() * i;
            count_num += partition_by_ucn[i].size();
        }
        else
        {
            count_num += (K - cur_uncon_num - cur_add_uc) / i;
            flag = true;
            break;
        }
    }
    if (count_nwoe != 0 && cur_add_uc < K && !flag)
    {
        int temp = 0;
        while (cur_add_uc < K && temp < count_nwoe)
        {
            cur_add_uc += count_num;
            count_num++;
            temp++;
        }
        if (cur_add_uc > K)
            count_num--;
    }
    if (min_degree_node.size() == 0 && min_degree == node_num)
    {
        min_degree_node.push_back(S_end);
    }
    assert(min_degree_node.size());
    ui rand_node = rand() % min_degree_node.size();
    swap_pos(S_end, min_degree_node[rand_node]);
    ui count_num2 = count_num;

    if (count_num > k_defective.size())
    {
        count_num2 = S_end;
        ui cur_k = K - cur_uncon_num;
        map<ui, ui> addk_num_pair;
        for (ui i = 0; i <= K; i++)
        {
            vector<ui> &nodes = partition_by_ucn[i];
            int color_num = cal_color_num(nodes);
            int remain = partition_by_ucn[i].size();
            ui addk = i;
            while (remain > 0)
            {
                if (addk_num_pair.find(addk) == addk_num_pair.end())
                {
                    addk_num_pair[addk] = min(color_num, remain);
                }
                else
                {
                    addk_num_pair[addk] += min(color_num, remain);
                }
                remain -= color_num;
                addk++;
            }
        }
        if (addk_num_pair.find(0) != addk_num_pair.end())
        {
            count_num2 += addk_num_pair[0];
        }
        for (ui i = 1; i <= K; i++)
        {
            if (addk_num_pair.find(i) == addk_num_pair.end())
                continue;
            if (addk_num_pair[i] * i >= cur_k)
            {
                return count_num2 + cur_k / i;
            }
            else
            {
                cur_k -= addk_num_pair[i] * i;
                count_num2 += addk_num_pair[i];
            }
        }
        if (count_nwoe != 0)
        {
            int temp = 0;
            while (cur_k >= count_num2 && temp < count_nwoe)
            {
                cur_k -= count_num2;
                count_num2++;
                temp++;
            }
        }
    }

    return min(count_num, count_num2);
}

void Graph::move_u_to_S(ui S_end, ui R_end, ui u)
{
    cur_uncon_num += S_end - degree_in_S[u];
    char *matrix_t = matrix + node_num * u;
    ui neighbors_n = 0;
    assert(SR_rid[u] == S_end);
    for (ui i = 0; i < R_end; i++)
    {
        if (i == S_end)
            continue;
        if (matrix_t[SR[i]])
        {
            degree_in_S[SR[i]]++;
            neighbors[neighbors_n++] = SR[i];
        }
    }
    for (ui i = 0; i < neighbors_n; i++)
    {
        for (ui j = i + 1; j < neighbors_n; j++)
        {
            assert(cn[neighbors[i] * node_num + neighbors[j]]);
            --cn[neighbors[i] * node_num + neighbors[j]];
            --cn[neighbors[j] * node_num + neighbors[i]];
        }
    }
}

void Graph::remove_u_from_S(ui S_end, ui R_end, ui u, ui level)
{
    swap_pos(S_end - 1, R_end - 1);
    levels[u] = level;
    char *matrix_t = matrix + node_num * u;
    for (ui i = 0; i < R_end - 1; i++)
    {
        if (matrix_t[SR[i]])
        {
            degree_in_S[SR[i]]--;
            degree[SR[i]]--;
        }
    }
    cur_uncon_num -= S_end - 1 - degree_in_S[u];
}

void Graph::restore(ui S_end, ui &R_end, ui level, ui old_r_end, ui old_edge_n)
{
    for (ui i = R_end; i < old_r_end; i++)
    {
        ui u = SR[i];
        assert(levels[u] == level + 1);
        ui neighbor_n = 0;
        char *matrix_t = matrix + u * node_num;
        degree[u] = 0;
        for (ui j = 0; j < R_end; j++)
            if (matrix_t[SR[j]])
            {
                ui w = SR[j];
                neighbors[neighbor_n++] = w;
                ++degree[w];
                ++degree[u];
            }
        for (ui j = 0; j < neighbor_n; j++)
        {
            ui v = neighbors[j];
            for (ui k = j + 1; k < neighbor_n; k++)
            {
                ui w = neighbors[k];
                ++cn[v * node_num + w];
                ++cn[w * node_num + v];
            }
        }
        ui *cn_t = cn + u * node_num;
        for (ui j = 0; j < R_end; j++)
            cn_t[SR[j]] = 0;
        for (ui j = 0; j < neighbor_n; j++)
            if (SR_rid[neighbors[j]] >= S_end)
            {
                ui v = neighbors[j];
                char *matrix_t = matrix + v * node_num;
                for (ui k = 0; k < R_end; k++)
                    if (matrix_t[SR[k]])
                        ++cn_t[SR[k]];
            }
        for (ui j = 0; j < R_end; j++)
        {
            cn[SR[j] * node_num + u] = cn_t[SR[j]];
        }
        ++R_end;
    }
    for (ui i = old_edge_n; i < removed_edge_n; i++)
    {
        ui v = removed_edges[i].first, w = removed_edges[i].second;
        assert(SR_rid[v] >= S_end && SR_rid[v] < R_end && SR_rid[w] >= S_end && SR_rid[w] < R_end);
        if (matrix[v * node_num + w])
            continue;
        matrix[v * node_num + w] = matrix[w * node_num + v] = 1;
        ++degree[v];
        ++degree[w];
        char *matrix_t = matrix + v * node_num;
        for (ui j = 0; j < R_end; j++)
            if (matrix_t[SR[j]])
            {
                ++cn[w * node_num + SR[j]];
                ++cn[SR[j] * node_num + w];
            }
        matrix_t = matrix + w * node_num;
        for (ui j = 0; j < R_end; j++)
            if (matrix_t[SR[j]])
            {
                ++cn[v * node_num + SR[j]];
                ++cn[SR[j] * node_num + v];
            }
    }
    removed_edge_n = old_edge_n;
}

void Graph::branch_and_bound(ui S_end, ui &R_end, ui v, ui level)
{
    branch_counting++;
    ui lb = k_defective.size();
    if (lb == node_num)
    {
        return;
    }
    if (cur_uncon_num > K)
    {
        return;
    }
    if (v != -1)
    {
        reduction(S_end, R_end, v, level);
    }
    if (R_end <= lb)
    {
        return;
    }
    if (S_end == R_end)
    {
        cout << "  find a better result of size:" << S_end << endl;
        k_defective.clear();
        for (ui i = 0; i < S_end; i++)
        {
            k_defective.push_back(node_trans[SR[i]]);
        }
        return;
    }
    ui ub;
    ub = cal_candidate_bound_colorset(S_end, R_end);
    ui old_r_end = R_end;
    ui old_edge_num = removed_edge_n;
    if (ub > lb)
    {
        ui v = SR[S_end];
        move_u_to_S(S_end, R_end, v);
        S_end++;
        branch_and_bound(S_end, R_end, v, level + 1);
        restore(S_end, R_end, level, old_r_end, old_edge_num);
        remove_u_from_S(S_end, R_end, v, level + 1);
        S_end--;
        R_end--;
        branch_and_bound(S_end, R_end, v, level + 1);
        restore(S_end, R_end, level, old_r_end, old_edge_num);
    }
    return;
}

void Graph::output_one_kdefective()
{
    FILE *fout = Utility::open_file("kdefective.txt", "w");
    fprintf(fout, "%lu\n", k_defective.size());
    sort(k_defective.begin(), k_defective.end());
    for (ui i = 0; i < k_defective.size(); i++)
        fprintf(fout, " %u", k_defective[i]);
    fprintf(fout, "\n");
    fclose(fout);
}

bool Graph::verify_result()
{
    ui uncon_count = 0;
    memset(exist, false, sizeof(bool) * n);
    for (ui i = 0; i < k_defective.size(); i++)
    {
        ui u = k_defective[i];
        for (ui j = pstart[u]; j < pend[u]; j++)
        {
            exist[edges[j]] = true;
        }
        for (ui j = i + 1; j < k_defective.size(); j++)
        {
            ui v = k_defective[j];
            if (!exist[v])
            {
                uncon_count++;
                if (uncon_count > K)
                {
                    cout << "WA!!! Unconnected node num exceed K" << endl;
                    return false;
                }
            }
        }
        for (ui j = pstart[u]; j < pend[u]; j++)
        {
            exist[edges[j]] = false;
        }
    }
    return true;
}

void Graph::k_defective_exact(FILE *f)
{
    cout << filename << endl;
    Timer t;
    k_defective.clear();
    fastLB2();
    if (k_defective.size() < 200 && t.elapsed() * 0.000001 < 200)
    {
        fastLB1();
    }
    ui R_end = preprocess();
    // output_pre_result();

    cout << "Branch and bound start!" << endl;
    Timer tt;
    branch_and_bound(0, R_end, -1, 0);
    searching_time = tt.elapsed() * 0.000001;
    cout << "Searching Time:" << searching_time << "s" << endl;
    result_num = k_defective.size();
    final_time = t.elapsed() * 0.000001;
    cout << "****** Total Time:" << final_time << "s"
         << " ******" << endl;
    fprintf(f, "%u\t%u\t%lf\t%lf\t%u\t%lf\n", pre_node, pre_edge, pre_time, searching_time, result_num, final_time);

    return;
}

void Graph::output_result()
{
    cout << "A k-defective of size " << k_defective.size() << " found : ";
    sort(k_defective.begin(), k_defective.end());
    for (ui i = 0; i < k_defective.size(); i++)
        cout << k_defective[i] << ' ';
    cout << endl;
    cout << "branching_times:" << branch_counting << endl;
}