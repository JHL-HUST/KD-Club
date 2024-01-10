#ifndef _GRAPH_H_
#define _GRAPH_H_

#include "Utility.h"
#include "Timer.h"
#include <queue>
#include <ctime>
#include <map>

class Graph
{
private:
    std::string filename; // input graph directory
    std::string fn;
    ui n;  // number of nodes of the graph
    ept m; // number of edges of the graph
    ui K;  // the value of k in k-plex
    std::vector<ui> k_defective;
    // For Preprocessing
    ept *pstart; // offset of neighbors of nodes
    ept *pend;   // used in search
    ept *pend2;
    ui *edges; // adjacent ids of edges
    ui *edges2;
    ui *degree;
    ui *tri_num;
    bool *exist;
    // For Branching
    char *matrix; // adjacent matrix
    ui *node_trans;
    ui *node_trans_rid;
    ui *cn; // common neighbors
    ui *degree_in_S;
    ui *SR;
    ui *SR_rid;
    ui *neighbors;
    ui *levels;
    ui node_num;
    ui cur_uncon_num;
    std::vector<std::pair<ui, ui>> removed_edges;
    std::vector<std::pair<ui, ui>> vp;
    ui removed_edge_n;
    ui rem_edge_num;

    double pre_time;
    ui pre_node;
    ui pre_edge;
    ui result_num;
    double final_time;
    double searching_time;

    ui branch_counting;

private:
    // heuristic algorithm for calculating LB
    int fastLB();
    ui degree_based_heu(ui *nodes, bool *in_cand, ui &s_end, ui &r_end);
    bool clique_assert();
    // algorithms for preprocessing
    ui check_vertex(ui *Qv, ui &Qv_n, ui lb);
    void check_edges(ui *Qe, ui &Qe_n, ui lb, ui rem_node_num);
    void triangle_count();
    void update_tri_num(ui u, ui v);
    ui construct_bnb_solver(ui rem_node_num);
    // algorithms for branching
    ui cal_color_num(std::vector<ui> &nodes);
    void reduction(ui S_end, ui &R_end, ui last_node, ui level);
    void move_u_to_S(ui S_end, ui R_end, ui u);
    void remove_u_from_S(ui S_end, ui R_end, ui u, ui level);
    void restore(ui S_end, ui &R_end, ui level, ui old_r_end, ui old_edge_n);
    ui cal_candidate_bound(ui S_end, ui R_end);
    ui cal_candidate_bound_colorset2(ui S_end, ui R_end);
    ui preprocess();
    void branch_and_bound(ui S_end, ui &R_end, ui v, ui level);
    void swap_pos(ui i, ui j);
    void show(ui S_end, ui R_end);

public:
    Graph(const char *_filename, const int _K);
    ~Graph();
    void read_graph();
    void k_defective_exact(FILE *f);
    void output_result();
    void output_one_kdefective();
    bool verify_result();

    // only for baseline
    void output_pre_result();
};
#endif