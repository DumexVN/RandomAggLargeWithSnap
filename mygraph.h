#ifndef GRAPH_H
#define GRAPH_H

#include <QtCore>
#include <QTimer>

#include "vertex.h"
#include "edge.h"

#include "Snap.h"

class Graph
{
public:
    Graph();
    //for working with snap
    PUNGraph convertToSnapUnGraph() const;
    bool convertCommToSnapVec(TVec < TCnCom > &CommV);
    bool convertSnapCommtoMyComm(const TVec < TCnCom > &CommV, QList<QList<quint32> > &result);
    //some simple Gnp generator
    void generateHiddenGnp(double pin, double pout);
    //
    void read_GML_file(QString filePath);
    void load_LFR_groundTruth();
    void read_Gephi_graph_and_produce_super_graph(QString filePath);
    void Gephi_parse_ModularityClass(int noClass);
    void save_edge_file_from_GML();
    void save_current_run_as_edge_file(QString fileName);
    void save_current_run_summary_file(QString fileName);
    void read_DUMEX_input(QString dirPath);
    void read_simple_edge(QString dirPath);
    void load_ground_truth_communities();
    void read_large_graph_with_ground_truth_communities(QString filePath);
    void read_superGraph(QString edgePath, QString summaryPath);
    void read_SuperCluster(QString filePath);
    //investigate bridges
    void get_bridge_stats();
    void LARGE_rerun();
    double LARGE_compute_modularity();
    double LARGE_compute_modularit_for_truth();
    //run
    void run_aggregation_on_selection(int n);
    void LARGE_reset();
    //stats
    double cal_average_clustering_coefficient();
    void clear_log();
    //Girvan and Newman Betweenness Centrality
    void betweenness_centrality_clustering();
    //aggregation
    void random_aggregate();
    void random_aggregate_with_degree_comparison();
    void random_aggregate_with_weight_comparison();
    void random_aggregate_with_neighbour_initial_degree_bias();
    void random_aggregate_with_neighbour_initial_degree_bias_with_constraint();
    void random_aggregate_with_neighbour_CURRENT_degree_bias();
    void random_aggregate_with_neighbour_CURRENT_degree_bias_with_constraint();
    void random_aggregate_highest_CURRENT_degree_neighbour();
    void random_aggregate_with_minimum_weight_neighbour();
    void random_aggregate_probabilistic_lowest_degree_neighbour_destructive();
    void random_aggregate_probabilistic_candidate_with_minimum_weight_neighbour();
    void random_aggregate_with_highest_edge_weight_and_weight_comparison();
    void random_aggregate_with_edge_weight_bias_and_weight_comparison();
    void random_aggregate_with_highest_triangulated_vertex();
    void random_aggregate_greedy_max_degree();
    void random_aggregate_greedy_max_weight();
    //agg without 'removing' vertices
    void random_aggregate_retain_vertex_using_triangulation();
    void random_aggregate_retain_vertex_using_probabilistic_triangulation();
    void random_aggregate_retain_vertex_using_triangulation_times_weight();
    void random_aggregate_retain_vertex_using_triangulation_of_cluster();
    //post aggregation
    void PostAgg_generate_super_vertex();
    void PostAgg_adjust_variables();
    void merge_result_clusters(const int &limit);
    void manual_set_working_dir(QString dirPath);
    void ReAgg_print_communities_stats();
    void ReAgg_select_and_save_community(int k);
    void ReAgg_select_and_save_community(const QList<quint32> &list);
    void ReAgg_select_and_save_community_with_intra_edge(const quint32 &size_minuimum);
    void ReAgg_select_and_save_community_by_size(const quint32 &size_minuimum);
    //restructing tree
    void construct_graph_at_a_specific_level(int level);

    //colouring
    void colouring_cluster_result();
    //quality
    QList<double> LARGE_compute_Pairwise_efficient(int n);
    double LARGE_compute_Newman_fraction_of_classified();
private:
    void assign_vertex_to_its_ground_truth_comm();
    bool locate_file_in_dir(QString &fileName);
    int get_number_from_filename(QString filename);
    void generate_base_graph_file(QString dirPath);
    QList<quint32> RGB_converter(quint32 hex);
    void read_ground_truth_communities();
    bool checkGraphCondition();
    void reConnectGraph();
    // for large graph
    void reindexing();
    void reindexing_ground_truth();
    void read_large_ground_truth_communities();
    void remove_excluded_vertices_from_ground_truth();
    void large_process_overlap();
    void large_process_overlap_by_seperate_intersection();
    void large_process_overlap_by_merge_intersection();
    void large_graph_parse_result();
    void large_parse_retain_result();
    void parse_LFR_groundTruth(QString filepath, int level);
    void print_result_stats();
    void print_single_community_inGraphML(const QList<quint32> &com, int k);
    void print_multiple_communities_inGraphML(const QList<quint32> &list);
    void print_multiple_communities_inGraphML(QString str);
    void print_single_community_with_only_intra_edge_inGraphML(int k);
    void print_result_community__with_attributes_inGraphML();
    void print_result_community__with_attributes_inGraphML(const QMap<quint32,QString> &colour);
    void LARGE_compute_cluster_matching(quint32 n);

    //colouring clusters
    void project_higher_levels_on_lower_levels(const QList<QPair<QString,QString> > &toParse);
    void mapping_colour_to_cluster(QMap<quint32, QString> &colourmap);


    bool LARGE_reload();
    void LARGE_reload_edges();
    void LARGE_reload_superEdges();
    void save_current_clusters();

    void parse_tree(const QList<QString> &file);
    void get_one_level_cluster(QString filename, QList<QList<quint32> > &clus);

    quint32 count_unique_element();
    quint64 calA(QList<quint64> param);
    quint64 calB(QList<quint64> param);
    quint64 calC(QList<quint64> param);
    quint64 calD(QList<quint64> param);
    double calAdRand(QList<quint64> param);

    //
    QString GMLpath;
    //
    QList<Vertex*> myVertexList;
    QList<Edge*> myEdgeList;
    QList<Vertex*> centroids;
    //
    QList<QList<quint32> > ground_truth_communities;
    QList<QPair<quint32,quint32> > hierarchy;
    QList<QList<quint32> > large_result;
    QSet<quint32> large_excluded;
    QMap<quint32, quint32> overlapped_vertices_ground_truth_cluster;
    //
    bool graphIsReady;
};
#endif // GRAPH_H
