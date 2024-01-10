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

    node_trans = nullptr;
    node_trans_rid = nullptr;
    matrix = nullptr;
    cn = nullptr;
    degree_in_S = nullptr;
    SR = nullptr;
    SR_rid = nullptr;
    neighbors = nullptr;
    levels = nullptr;
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
    degree = new ui[n];
    // pstart记录每一个点在vp中开始的位置（主要记录每个点在edges的范围）
    // edges记录当前vp中边的另一个端点
    pstart[0] = 0;
    ui idx = 0;
    for (ui i = 0; i < n; i++)
    {
        pstart[i + 1] = pstart[i];
        while (idx < vp.size() && vp[idx].first == i)
            edges[pstart[i + 1]++] = vp[idx++].second;
        degree[i] = pstart[i + 1] - pstart[i];
    }

    fclose(f);

#ifndef NDEBUG
    printf("Finished reading graph\n");
#endif
}

int Graph::fastLB()
{
    cout << "Start FastLB to calculate a lower bound!" << endl;
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
        ui cli_size = degree_based_heu(nodes, in_cand, s_end, r_end);
        if (cli_size > max_clique)
        {
            k_defective.clear();
            for (ui i = 0; i < s_end; i++)
            {
                k_defective.push_back(nodes[i]);
            }
            max_clique = cli_size;
        }
    }
    if (!clique_assert())
    {
        cout << "verify clique failed, please check the codes!" << endl;
    }
    cout << "FastLB finishes! Find a clique of size " << k_defective.size() << " at time " << Utility::integer_to_string(t.elapsed()).c_str() << "." << endl;
    delete[] nodes;
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
    return degree_based_heu(nodes, in_cand, s_end, r_end);
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
            if (v > u)
            {
                edges2[pend2[v]] = u;
                tri_num[pend2[v]] = tri_num[i];
                pend2[v]++;
            }
        }
    }
}
ui Graph::check_vertex(ui *Qv, ui &Qv_n, ui lb)
{
    cout << "  Start checking useless nodes!" << endl;
    int rem_node_num = n;
    ui removed_node_count = 0;
    bool *in_queue = new bool[n];
    memset(in_queue, false, sizeof(bool) * n);
    for (ui u = 0; u < n; u++)
    {
        ui ub = 1 + degree[u] + min(rem_node_num - degree[u] - 1, K);
        if (ub <= lb)
        {
            exist[u] = false;
            removed_node_count++;
            rem_node_num--;
            for (ui j = pstart[u]; j < pstart[u + 1]; j++)
            {
                ui v = edges[j];
                if (!exist[v])
                {
                    continue;
                }
                degree[v]--;
                ui ub_v = 1 + degree[v] + min(rem_node_num - degree[v] - 1, K);
                if (v < u && !in_queue[v] && ub_v <= lb)
                {
                    Qv[Qv_n++] = v;
                    in_queue[v] = true;
                }
            }
        }
    }
    for (ui i = 0; i < Qv_n; i++)
    {
        ui u = Qv[i];
        assert(in_queue[u] && exist[u]);
        exist[u] = false;
        removed_node_count++;
        rem_node_num--;
        for (ui j = pstart[u]; j < pstart[u + 1]; j++)
        {
            ui v = edges[j];
            if (!exist[v])
            {
                continue;
            }
            degree[v]--;
            ui ub_v = 1 + degree[v] + min(rem_node_num - degree[v] - 1, K);
            if (!in_queue[v] && ub_v <= lb)
            {
                Qv[Qv_n++] = v;
                in_queue[v] = true;
            }
        }
        in_queue[u] = false;
    }
    for (ui u = 0; u < n; u++)
    {
        pend2[u] = pstart[u];
        pend[u] = pstart[u + 1];
        if (!exist[u])
        {
            degree[u] = 0;
            continue;
        }
        for (ui j = pstart[u]; j < pend[u]; j++)
        {
            ui v = edges[j];
            if (exist[v] && u < v)
            {
                edges2[pend2[u]++] = v;
            }
        }
    }
    delete[] in_queue;
    cout << "    After checking vertex, " << removed_node_count << " nodes are removed! " << rem_node_num << " nodes remain!" << endl;
    return rem_node_num;
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
void Graph::check_edges(ui *Qe, ui &Qe_n, ui lb, ui rem_node_num)
{
    cout << "  Start checking useless edges!" << endl;
    triangle_count();
    ui reomved_edge_count = 0;
    bool *in_queue = new bool[n];
    memset(in_queue, false, sizeof(bool) * n);
    for (ui u = 0; u < n; u++)
    {
        if (!exist[u])
            continue;
        for (ui j = pstart[u]; j < pend2[u]; j++)
        {
            ui v = edges2[j];
            if (!exist[v] || v < u)
                continue;
            ui ub = 2 + tri_num[j] + min(K, rem_node_num - tri_num[j] - 2);
            if (ub <= lb)
            {
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

                if (!in_queue[v])
                {
                    Qe[Qe_n++] = v;
                    in_queue[v] = true;
                }
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
            if (!exist[v])
                continue;

            ui ub = 2 + tri_num[j] + min(K, rem_node_num - tri_num[j] - 2);
            if (ub <= lb)
            {
                remove = true;
                --pend2[u];
                edges2[j] = edges2[pend2[u]];
                tri_num[j] = tri_num[pend2[u]];
                --j;
                reomved_edge_count++;
                update_tri_num(u, v);
                if (!in_queue[v])
                {
                    Qe[Qe_n++] = v;
                    in_queue[v] = true;
                }
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
    memset(exist, true, sizeof(bool) * n);
    ui Qv_n = 0;
    ui *Qv = new ui[n];
    ui lb = k_defective.size();
    ui rem_node_num = check_vertex(Qv, Qv_n, lb);
    delete[] Qv;
    ui Qe_n = 0;
    ui *Qe = new ui[m];
    check_edges(Qe, Qe_n, lb, rem_node_num);
    delete[] Qe;
    memset(degree, 0, sizeof(ui) * n);
    for (ui u = 0; u < n; u++)
    {
        if (exist[u])
        {
            degree[u] = pend2[u] - pstart[u];
        }
    }
    for (ui u = 0; u < n; u++)
    {
        if (!exist[u])
            continue;
        ui ub = 1 + degree[u] + min(rem_node_num - degree[u] - 1, K);
        if (ub <= lb)
        {
            exist[u] = false;
            rem_node_num--;
            for (ui j = pstart[u]; j < pend2[u]; j++)
            {
                ui v = edges2[j];
                if (!exist[v])
                {
                    continue;
                }
                degree[v]--;
            }
        }
    }
    rem_edge_num = 0;
    for (ui u = 0; u < n; u++)
    {
        if (!exist[u])
            continue;
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

    string wf = "./file/" + fn.substr(0, fn.find('.')) + "-" + to_string(K) + ".col";
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
    memset(degree_in_S, 0, sizeof(ui) * rem_node_num);
    memset(matrix, 0, sizeof(char) * rem_node_num * rem_node_num);
    memset(cn, 0, sizeof(ui) * rem_node_num * rem_node_num);
    memset(degree, 0, sizeof(ui) * rem_node_num);
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
            R_end++;
        }
    }
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
        ui v = SR[i];
        if (matrix_t[v])
        {
            Q[Q_n++] = v;
            in_queue[v] = true;
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
    for (ui i = 0; i < Q_n; i++)
    {
        ui w = Q[i];
        if (SR_rid[w] >= R_end)
        {
            continue;
        }
        ui *cn_t = cn + w * node_num;
        matrix_t = matrix + w * node_num;
        bool on_show = false;
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
                matrix[v * node_num + w] = 0;
                matrix[w * node_num + v] = 0;
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

ui Graph::cal_candidate_bound(ui S_end, ui R_end)
{
    int counter[K + 1];
    for (int i = 0; i <= K; i++)
    {
        counter[i] = 0;
    }
    ui min_degree = node_num;
    vector<ui> min_degree_node;
    for (ui i = S_end; i < R_end; i++)
    {
        int uc_num = S_end - degree_in_S[SR[i]];
        assert(uc_num <= K && uc_num >= 0);
        counter[uc_num]++;
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
    for (int i = 0; i <= K; i++)
    {
        if (cur_add_uc + counter[i] * i <= K - cur_uncon_num)
        {
            cur_add_uc += counter[i] * i;
            count_num += counter[i];
        }
        else
        {
            count_num += (K - cur_uncon_num - cur_add_uc) / i;
            break;
        }
    }
    assert(min_degree_node.size());
    ui rand_node = rand() % min_degree_node.size();
    swap_pos(S_end, min_degree_node[rand_node]);
    return count_num;
}
ui Graph::cal_candidate_bound_colorset2(ui S_end, ui R_end)
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
    ui ub = cal_candidate_bound(S_end, R_end);
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
    ui lb = fastLB();
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