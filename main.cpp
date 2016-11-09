#include <iostream>

#include <QCoreApplication>
#include <QtGlobal>
#include <QtDebug>
#include <QFile>
#include <QTextStream>

#include <random>
#include <stdio.h>
#include <stdlib.h>
#include "mygraph.h"

#include "Snap.h"

#include "boost/graph/connected_components.hpp"
#include <boost/graph/adjacency_list.hpp>

QString workingDir = "C:/Users/Dumex/Desktop/"; //change this

enum RandomAgg { I_a, I_b, I_c,
                 II_a, II_a_i, II_b, II_b_i, II_c, II_d, II_e, II_f, II_g, II_h,
                 III_a, III_b, III_c, III_d, III_e,
                 GN_Clustering, CNM_Clustering,
                 R1a,
                 I_x,
                 RFD,
                 III_z
               };
const char *name[] = { "I.a", "I.b", "I_c",
                       "II.a", "II.a.i", "II.b", "II.b.i", "II.c", "II.d", "II.e", "II.f", "II.g", "II.h",
                       "III.a", "III.b", "III.c"," III.d", "III.e",
                       "GN_Clustering", "CNM_Clustering",
                       "R1a",
                       "I_x",
                       "RFD",
                       "III.z"
                     };
/** Override Debug Message Handler
 * @brief myMessageOutput
 * @param type
 * @param context
 * @param msg
 */
void myMessageOutput(QtMsgType type, const QMessageLogContext &context, const QString &msg)
{
    QString txt;
    switch (type) {
    case QtDebugMsg:
        txt =  QString("%1 ").arg(msg).append('/n');
        break;
    case QtFatalMsg:
        txt =  QString("Fatal: %1 ").arg(msg).append('/n');
        abort();
    }
    QFile outFile("C:/Users/Dumex/Desktop/SocialNetworksCollection/cond-mat-2003/log.txt");
    outFile.open(QIODevice::WriteOnly | QIODevice::Append);
    QTextStream ts(&outFile);
    ts << txt << endl;
    outFile.close();
}

/** Overload Since Rebuilding Graph Requires Re-read the edgefile.txt
 * @brief writeEdgeFile
 */
void writeEdgeFile(const Graph &G)
{
    PUNGraph TGraph = G.convertToSnapUnGraph();
    QString filePath = workingDir + "edge_file.txt";
    std::string dir = filePath.toStdString();
    char * c = &dir[0];
    TStr str(c);
    TSnap::SaveEdgeList(TGraph,
                        str,
                        "Save as tab-separated list of edges");
}

void writeSeperateFile(QString fileName, const QList<QList<double> > &matrix)
{
    qDebug() << "Writing File";
    QFile file(workingDir + fileName + ".txt");
    file.open(QFile::WriteOnly | QFile::Text);
    QTextStream out(&file);
    out.setRealNumberPrecision(4);
    int h = matrix.size(), w = matrix[0].size();
    for (int i = 0; i < h; i++)
    {
        QString str;
        for (int j = 0; j < w; j++)
        {
            str.append(QString::number(matrix[i][j])).append('/t');
        }
        out << str;
        out << '/n';
    }
    file.close();
}

void plotMatrix(const QList<QList<double> > &matrix, QString str)
{
    TVec<TVec<TPair<TInt, TFlt > > > AllSeries;
    int h = matrix.size(), w = matrix[0].size();
    std::string stdstr = str.toStdString();
    char * c = &stdstr[0];
    TStr label(c);
    //constructing each series
    for (int j = 0; j < w; j++)
    {
        TVec<TPair<TInt, TFlt> > XY;
        for(int i = 0; i < h; i++)
        {
            TPair<TInt, TFlt> p(i, matrix[i][j]); //for continous x i.e. GN Experiment
          //  TPair<TInt, TFlt> p(i+5, matrix[i][j]);
            XY.Add(p);
        }
        AllSeries.Add(XY);
    }
    //set plot desription
    TStr plot_desc = c;
    plot_desc+= "plot";
    TVec<TGnuPlot> plotV;
    for(int i = 0; i < 4; i++)
    {
        TStr lab = label;
        if (i == 0)lab+= "_TypeI";
        else if (i == 1)lab+= "_TypeII";
        else if (i == 2)lab+= "_TypeIII";
        else lab+= "_Best";
        TGnuPlot Gp(lab, plot_desc);
        Gp.SetYRange(-0.5,1.0);
        plotV.Add(Gp);
    }
    //plot
    for(int i = 0; i < AllSeries.Len(); i++)
    {
        if (i <= RandomAgg::I_c || i == RandomAgg::R1a || i == RandomAgg::I_x) plotV[0].AddPlot(AllSeries[i], gpwLinesPoints, name[i]);
        else if (i > RandomAgg::I_c && i <= RandomAgg::II_h) plotV[1].AddPlot(AllSeries[i], gpwLinesPoints, name[i]);
        else if (i > RandomAgg::II_h && i <= RandomAgg::III_e) plotV[2].AddPlot(AllSeries[i], gpwLinesPoints, name[i]);
        else {}
    }
    //plot the best
    plotV[3].AddPlot(AllSeries[RandomAgg::I_c], gpwLinesPoints, name[RandomAgg::I_c]);
    plotV[3].AddPlot(AllSeries[RandomAgg::II_g], gpwLinesPoints, name[RandomAgg::II_g]);
    plotV[3].AddPlot(AllSeries[RandomAgg::III_a], gpwLinesPoints, name[RandomAgg::III_a]);
    plotV[3].AddPlot(AllSeries[RandomAgg::III_e], gpwLinesPoints, name[RandomAgg::III_e]);
    plotV[3].AddPlot(AllSeries[RandomAgg::GN_Clustering], gpwLinesPoints, name[RandomAgg::GN_Clustering]);
    plotV[3].AddPlot(AllSeries[RandomAgg::CNM_Clustering], gpwLinesPoints, name[RandomAgg::CNM_Clustering]);
    plotV[3].AddPlot(AllSeries[RandomAgg::RFD], gpwLinesPoints, name[RandomAgg::RFD]);
    plotV[3].AddPlot(AllSeries[RandomAgg::III_z], gpwLinesPoints, name[RandomAgg::III_z]);

    //export as png
    for (int i = 0; i < plotV.Len();i++)
    {
        TGnuPlot Gp = plotV[i];
        Gp.SetXYLabel("Z_out", label);
        Gp.SavePng(label+".png");
    }
}


void plotIntMatrix(const QList<QList<int> > &matrix, QString str)
{
    TVec<TVec<TPair<TInt, TInt > > > AllSeries;
    int h = matrix.size(), w = matrix[0].size();
    std::string stdstr = str.toStdString();
    char * c = &stdstr[0];
    TStr label(c);
    //constructing each series
    for (int j = 0; j < w; j++)
    {
        TVec<TPair<TInt, TInt> > XY;
        for(int i = 0; i < h; i++)
        {
            TPair<TInt, TInt> p((i+1)*5, matrix[i][j]); //for continous x i.e. GN Experiment
            XY.Add(p);
        }
        // printf("XY:%d/n", XY.Len() );
        AllSeries.Add(XY);
    }
    //set plot desription
    TStr plot_desc = c;
    plot_desc+= "plot";
    TVec<TGnuPlot> plotV;
    for(int i = 0; i < 4; i++)
    {
        TStr lab = label;
        if (i == 0)lab+= "_TypeI";
        else if (i == 1)lab+= "_TypeII";
        else if (i == 2)lab+= "_TypeIII";
        else lab+= "_Best";
        TGnuPlot Gp(lab, plot_desc);
        plotV.Add(Gp);
    }
    //plot
    plotV[3].AddPlot(AllSeries[0], gpwLinesPoints, name[RandomAgg::RFD]);
    /*
    for(int i = 0; i < AllSeries.Len(); i++)
    {
        if (i <= RandomAgg::I_c) plotV[0].AddPlot(AllSeries[i], gpwLinesPoints, name[i]);
        else if (i > RandomAgg::I_c && i <= RandomAgg::II_h) plotV[1].AddPlot(AllSeries[i], gpwLinesPoints, name[i]);
        else if (i > RandomAgg::II_h && i <= RandomAgg::III_e) plotV[2].AddPlot(AllSeries[i], gpwLinesPoints, name[i]);
        else
        {
          //  plotV[3].AddPlot(AllSeries[RandomAgg::I_c], gpwLinesPoints, name[RandomAgg::I_c]);
            plotV[3].AddPlot(AllSeries[RandomAgg::RFD], gpwLinesPoints, name[RandomAgg::RFD]);
         //   plotV[3].AddPlot(AllSeries[RandomAgg::II_g], gpwLinesPoints, name[RandomAgg::II_g]);
         //   plotV[3].AddPlot(AllSeries[RandomAgg::III_a], gpwLinesPoints, name[RandomAgg::III_a]);
         //   plotV[3].AddPlot(AllSeries[RandomAgg::III_e], gpwLinesPoints, name[RandomAgg::III_e]);
         //   plotV[3].AddPlot(AllSeries[RandomAgg::GN_Clustering], gpwLinesPoints, name[RandomAgg::GN_Clustering]);
         //   plotV[3].AddPlot(AllSeries[RandomAgg::CNM_Clustering], gpwLinesPoints, name[RandomAgg::CNM_Clustering]);
        }
    }
    //plot the best*/

    //export as png
    for (int i = 0; i < plotV.Len();i++)
    {
        TGnuPlot Gp = plotV[i];
        Gp.SetXYLabel("Out_Degree", label);
        Gp.SavePng();
    }
}

/** Girvan and Newman Experiments on 4 hidden partition of Gnp
 * @brief GN_experiment
 */
void GN_experiment()
{
    int n = 128, m = 4, n_per_c = n/m, neighbours = n_per_c - 1, outer = n_per_c*3;
    QList<QList<double> > RAND, JACCARD, ARI, Q, GN;
    for (int z_out = 0 ; z_out <= 12; z_out++)
    {
        QList<double> sRAND, sJACCARD, sARI, sQ, sGN; //s = Set
        for (int k = 0; k <= 22; k++)
        {
            //generate GN graph
            //run algorithm here
            double p_in = 0.0, p_out = 0.0;
            p_out = (double) z_out/outer;
            p_in = (double) (16-z_out)/neighbours;
            int times = 50;
            double iRAND = 0.0 , iJACCARD = 0.0, iARI = 0.0, iQ = 0.0, iGN = 0.0; // i = iterator
            QString mess;
            if (k == 0) {mess.append(QString("********** I.a ************ /n"));}
            else if (k == 1) {mess.append(QString( "********** I.b ************ /n"));}
            else if (k == 2) {mess.append(QString( "********** I.c ************ /n"));}
            else if (k == 3) {mess.append(QString( "********** II.a ************ /n"));}
            else if (k == 4){mess.append(QString( "********** II.a(i) ************ /n"));}
            else if (k == 5) {mess.append(QString( "********** II.b ************ /n"));}
            else if (k == 6){mess.append(QString( "********** II.b(i) ************ /n"));}
            else if (k == 7) {mess.append(QString( "********** II.c ************ /n"));}
            else if (k == 8) {mess.append(QString( "********** II.d ************ /n"));}
            else if (k == 9) {mess.append(QString( "********** II.e ************ /n"));}
            else if (k == 10) {mess.append(QString( "********** II.f ************ /n"));}
            else if (k == 11) {mess.append(QString( "********** II.h ************ /n"));}
            else if (k == 12){mess.append(QString( "********** II.g ************ /n"));}
            else if (k == 13){mess.append(QString( "********** III.a ************ /n"));}
            else if (k == 14){mess.append(QString( "********** III.b ************ /n"));}
            else if (k == 15){mess.append(QString( "********** III.c ************ /n"));}
            else if (k == 16){mess.append(QString( "********** III.d ************ /n"));}
            else if (k == 17){mess.append(QString( "********** III.e ************ /n"));}
            else if (k == 18){mess.append(QString( "********** Betweenness Centrality Clustering ************ /n")); times = 1;}
            else if (k == 19){mess.append(QString( "********** CNM Clustering ************ /n")); times = 1;}
            else if (k == 20){mess.append(QString( "********** R1a ************ /n"));}
            else if (k == 21){mess.append(QString( "********** I.x ************ /n"));}
            else if (k == 22){mess.append(QString( "********** RFD ************ /n"));}
            else{}
            for (int i = 0 ; i < times; i++)
            {
                Graph G;
                G.manual_set_working_dir(workingDir);
                G.generateHiddenGnp(p_in, p_out);
                writeEdgeFile(G);
                if (k == 0) {G.random_aggregate();}
                else if (k == 1) {G.random_aggregate_with_degree_comparison(); }
                else if (k == 2) {G.random_aggregate_with_weight_comparison(); }
                else if (k == 3) {G.random_aggregate_with_neighbour_initial_degree_bias(); }
                else if (k == 4){G.random_aggregate_with_neighbour_initial_degree_bias_with_constraint();}
                else if (k == 5) {G.random_aggregate_with_neighbour_CURRENT_degree_bias(); }
                else if (k == 6){G.random_aggregate_with_neighbour_CURRENT_degree_bias_with_constraint();}
                else if (k == 7) {G.random_aggregate_highest_CURRENT_degree_neighbour();}
                else if (k == 8) {G.random_aggregate_with_minimum_weight_neighbour();}
                else if (k == 9) {G.random_aggregate_probabilistic_lowest_degree_neighbour_destructive();}
                else if (k == 10) {G.random_aggregate_probabilistic_candidate_with_minimum_weight_neighbour();}
                else if (k == 11) {G.random_aggregate_greedy_max_weight();}
                else if (k == 12){G.random_aggregate_greedy_max_degree();}
                else if (k == 13){G.random_aggregate_retain_vertex_using_triangulation();}
                else if (k == 14){G.random_aggregate_retain_vertex_using_probabilistic_triangulation();}
                else if (k == 15){G.random_aggregate_with_highest_triangulated_vertex();}
                else if (k == 16){G.random_aggregate_retain_vertex_using_triangulation_times_weight();}
                else if (k == 17){G.random_aggregate_retain_vertex_using_triangulation_of_cluster();}
                else if (k == 18){G.betweenness_centrality_clustering();}
                else if (k == 19){G.fast_CMN();}
                else if (k == 20){G.reverse_random_aggregate();}
                else if (k == 21){G.reverse_random_aggregate_with_degree_comparison();}
                else if (k == 22){G.random_functional_digraph();}
                QList<double> id = G.LARGE_compute_Pairwise_efficient(-1); //param -1 set n to |V|
                iRAND+=id[0];
                iJACCARD+=id[1];
                iARI+=id[2];
                iQ+=G.LARGE_compute_modularity();
                double gn = G.compute_GN_index();
                G.LARGE_hard_reset();
                iGN += gn;
            }
            //add to list
            double normalisedRand = iRAND/times,
                    normalisedJaccard = iJACCARD/times,
                    normalisedARI = iARI/times,
                    normalisedQ = iQ/times,
                    normalisedGN = (double) iGN/times;
            sRAND << normalisedRand; sJACCARD << normalisedJaccard; sARI << normalisedARI; sQ << normalisedQ; sGN << normalisedGN;
            //write to file
            QFile file(workingDir +"Stat.txt");
            file.open(QFile::Text | QFile::Append);
            QTextStream out(&file);
            mess.append(QString("%1/t%2/t%3/t%4/t%5/t%6/n").arg(normalisedRand).arg(normalisedJaccard).arg(normalisedARI).arg(normalisedQ).arg(normalisedGN).arg(z_out));
            out << mess;
            file.close();
        }
        RAND << sRAND; JACCARD << sJACCARD; ARI << sARI; Q << sQ; GN << sGN;
    }
    writeSeperateFile(QString("ARI"), ARI);
    writeSeperateFile(QString("JACCARD"), JACCARD);
    writeSeperateFile(QString("Q"), Q);
    writeSeperateFile(QString("GirvanNewmanIndex"), GN);
    plotMatrix(ARI, QString("ARI"));
    plotMatrix(JACCARD, QString("Jaccard"));
    plotMatrix(Q, QString("Q"));
    plotMatrix(GN, QString("GirvanNewman"));
}

/** LARGE gnp experiment
 * See note for more details
 * q >> p
 * p < Tria_threshold derived from first principal (
 * (There are two thresholds for triangulation on Gnp:
 * 1. Vu's definition: on edge
 * 2. Colin's definition: on vertex
 * )
 * @brief LARGE_GN_experiment
 */
void LARGE_Gnp_experiment()
{
    int step = 10;
    quint32 ell = 1000, n = ell*4;
    double global_p_threshold = std::pow((double)n,(double)1/2);
    double layer_q = 0.6,
           p_stepsize = global_p_threshold/step;
    /*kout_threshold = std::pow((double)(ell*3),(double)1/2);
    double p_threshold = kout_threshold/(ell*3),
           q_threshold = kin_threshold/ell;
    double p_stepsize = p_threshold/step;*/


    //test_q = 0.6;

    QList<int> blocked_type;
//    blocked_type << 18 << 19 << 0 << 1 << 2 << 3 << 4 << 5 << 6 << 7 << 8 << 9 << 10 << 11 << 12;

    QList<QList<double> > RAND, JACCARD, ARI, Q, GN;
    for (double p = 0.0; p <= global_p_threshold; p += p_stepsize)
    {
        QList<double> sRAND, sJACCARD, sARI, sQ, sGN; //s = Set
        for (int k = 0; k <= 23; k++)
        {
            qDebug() << "-------------------------- New Algorithm -------------------- [p: " << p << " - thres:" << p_threshold << "]" ;
            //generate GN graph
            //run algorithm here
            int times = 3;
            double iRAND = 0.0 , iJACCARD = 0.0, iARI = 0.0, iQ = 0.0, iGN = 0.0; // i = iterator
            QString mess;
            if (k == 0) {mess.append(QString("********** I.a ************ \n"));}
            else if (k == 1) {mess.append(QString( "********** I.b ************ \n"));}
            else if (k == 2) {mess.append(QString( "********** I.c ************ \n"));}
            else if (k == 3) {mess.append(QString( "********** II.a ************ \n"));}
            else if (k == 4){mess.append(QString( "********** II.a(i) ************ \n"));}
            else if (k == 5) {mess.append(QString( "********** II.b ************ \n"));}
            else if (k == 6){mess.append(QString( "********** II.b(i) ************ \n"));}
            else if (k == 7) {mess.append(QString( "********** II.c ************ \n"));}
            else if (k == 8) {mess.append(QString( "********** II.d ************ \n"));}
            else if (k == 9) {mess.append(QString( "********** II.e ************ \n"));}
            else if (k == 10) {mess.append(QString( "********** II.f ************ \n"));}
            else if (k == 11) {mess.append(QString( "********** II.h ************ \n"));}
            else if (k == 12){mess.append(QString( "********** II.g ************ \n"));}
            else if (k == 13){mess.append(QString( "********** III.a ************ \n"));}
            else if (k == 14){mess.append(QString( "********** III.b ************ \n"));}
            else if (k == 15){mess.append(QString( "********** III.c ************ \n"));}
            else if (k == 16){mess.append(QString( "********** III.d ************ \n"));}
            else if (k == 17){mess.append(QString( "********** III.e ************ \n"));}
            else if (k == 18){mess.append(QString( "********** Betweenness Centrality Clustering ************ \n")); times = 1;}
            else if (k == 19){mess.append(QString( "********** CNM Clustering ************ \n")); times = 1;}
            else if (k == 20){mess.append(QString( "********** R1a ************ \n"));}
            else if (k == 21){mess.append(QString( "********** I.x ************ \n"));}
            else if (k == 22){mess.append(QString( "********** RFD ************ \n"));}
            else if (k == 23){mess.append(QString( "********** III.z ************ \n"));}
            else{}
            for (int i = 0 ; i < times; i++)
            {
                if (blocked_type.contains(k))
                {
                    iRAND+=0;
                    iJACCARD+=0;
                    iARI+=0;
                    iQ+=0;
                    iGN+=0;
                }
                else
                {
                    qDebug() << "----------- New Gnp ---------- [k: " << (k) << " ]";
                    Graph G;
                    G.manual_set_working_dir(workingDir);
                    G.generateHiddenGnp_LargeN_layered(p, layer_q, l);
                    writeEdgeFile(G);
                    if (k == 0) {G.random_aggregate();}
                    else if (k == 1) {G.random_aggregate_with_degree_comparison(); }
                    else if (k == 2) {G.random_aggregate_with_weight_comparison(); }
                    else if (k == 3) {G.random_aggregate_with_neighbour_initial_degree_bias(); }
                    else if (k == 4){G.random_aggregate_with_neighbour_initial_degree_bias_with_constraint();}
                    else if (k == 5) {G.random_aggregate_with_neighbour_CURRENT_degree_bias(); }
                    else if (k == 6){G.random_aggregate_with_neighbour_CURRENT_degree_bias_with_constraint();}
                    else if (k == 7) {G.random_aggregate_highest_CURRENT_degree_neighbour();}
                    else if (k == 8) {G.random_aggregate_with_minimum_weight_neighbour();}
                    else if (k == 9) {G.random_aggregate_probabilistic_lowest_degree_neighbour_destructive();}
                    else if (k == 10) {G.random_aggregate_probabilistic_candidate_with_minimum_weight_neighbour();}
                    else if (k == 11) {G.random_aggregate_greedy_max_weight();}
                    else if (k == 12){G.random_aggregate_greedy_max_degree();}
                    else if (k == 13){G.random_aggregate_retain_vertex_using_triangulation();}
                    else if (k == 14){G.random_aggregate_retain_vertex_using_probabilistic_triangulation();}
                    else if (k == 15){G.random_aggregate_with_highest_triangulated_vertex();}
                    else if (k == 16){G.random_aggregate_retain_vertex_using_triangulation_times_weight();}
                    else if (k == 17){G.random_aggregate_retain_vertex_using_triangulation_of_cluster();}
                    else if (k == 18){G.betweenness_centrality_clustering();}
                    else if (k == 19){G.fast_CMN();}
                    else if (k == 20){G.reverse_random_aggregate();}
                    else if (k == 21){G.reverse_random_aggregate_with_degree_comparison();}
                    else if (k == 22){G.random_functional_digraph();}
                    else if (k == 23){G.random_aggregate_retain_vertex_using_colin_triangulation();}
                    QList<double> id = G.LARGE_compute_Pairwise_efficient(n); //param -1 set n to |V|
                    iRAND+=id[0];
                    iJACCARD+=id[1];
                    iARI+=id[2];
                    iQ+=G.LARGE_compute_modularity();
                    double gn = 0.0;
    //                double gn = G.compute_GN_index();
                    G.LARGE_hard_reset();
                    iGN += gn;
                }
            }
            //add to list
            double normalisedRand = iRAND/times,
                    normalisedJaccard = iJACCARD/times,
                    normalisedARI = iARI/times,
                    normalisedQ = iQ/times,
                    normalisedGN = (double) iGN/times;
            sRAND << normalisedRand; sJACCARD << normalisedJaccard; sARI << normalisedARI; sQ << normalisedQ; sGN << normalisedGN;
            //write to file
            QFile file(workingDir +"Stat.txt");
            file.open(QFile::Text | QFile::Append);
            QTextStream out(&file);
            mess.append(QString("%1\t%2\t%3\t%4\t%5\t%6\n").arg(normalisedRand).arg(normalisedJaccard).arg(normalisedARI).arg(normalisedQ).arg(normalisedGN).arg(p));
            out << mess;
            file.close();
        }
        RAND << sRAND; JACCARD << sJACCARD; ARI << sARI; Q << sQ; GN << sGN;
    }
    writeSeperateFile(QString("ARI"), ARI);
    writeSeperateFile(QString("JACCARD"), JACCARD);
    writeSeperateFile(QString("Q"), Q);
    writeSeperateFile(QString("GirvanNewmanIndex"), GN);
    plotMatrix(ARI, QString("ARI"));
    plotMatrix(JACCARD, QString("Jaccard"));
    plotMatrix(Q, QString("Q"));
    plotMatrix(GN, QString("GirvanNewman"));
}

/////////////////////////////////////////////////////////////////////////////////////////
/// /brief Experiment with COMPLETE GRAPH
/// /param
///
void generateKN(int n)
{
    PUNGraph TGraph = PUNGraph::New();
    TGraph = TSnap::GenFull<PUNGraph>(n);
    TSnap::SaveEdgeList(TGraph,
                        "C:/Users/Dumex/Desktop/SocialNetworksCollection/CompleteGraph/edge_file.txt",
                        "K_n");
}

void CompleteGraph_Exp()
{
    Graph G;
    QFile file("C:/Users/Dumex/Desktop/SocialNetworksCollection/CompleteGraph/Stat.txt");
    QList<QList<double> > Q;
    QList<QList<int> > CC;
    for (int n = 1000; n <= 5000; n+=1000)
    {
        QList<double> sQ; //s = Set
        QList<int> sCC;
        //generate Watts Schorgat
        generateKN(n);
        G.LARGE_hard_reset();
        G.read_edge(QString("C:/Users/Dumex/Desktop/SocialNetworksCollection/CompleteGraph/"));
        for (int k = 0; k <= 19; k++)
        {
            //generate GN graph
            //run algorithm here
            int times = 10;
            double iQ = 0.0; // i = iterator
            quint32 iCC = 0;
            QString mess;
            if (k == 0) {mess.append(QString("********** I.a ************ /n"));}
            else if (k == 1) {mess.append(QString( "********** I.b ************ /n"));}
            else if (k == 2) {mess.append(QString( "********** I.c ************ /n"));}
            else if (k == 3) {mess.append(QString( "********** II.a ************ /n"));}
            else if (k == 4){mess.append(QString( "********** II.a(i) ************ /n"));}
            else if (k == 5) {mess.append(QString( "********** II.b ************ /n"));}
            else if (k == 6){mess.append(QString( "********** II.b(i) ************ /n"));}
            else if (k == 7) {mess.append(QString( "********** II.c ************ /n"));}
            else if (k == 8) {mess.append(QString( "********** II.d ************ /n"));}
            else if (k == 9) {mess.append(QString( "********** II.e ************ /n"));}
            else if (k == 10) {mess.append(QString( "********** II.f ************ /n"));}
            else if (k == 11) {mess.append(QString( "********** II.g ************ /n"));}
            else if (k == 12){mess.append(QString( "********** II.h ************ /n"));}
            else if (k == 13){mess.append(QString( "********** III.a ************ /n"));}
            else if (k == 14){mess.append(QString( "********** III.b ************ /n"));}
            else if (k == 15){mess.append(QString( "********** III.c ************ /n"));}
            else if (k == 16){mess.append(QString( "********** III.d ************ /n"));}
            else if (k == 17){mess.append(QString( "********** III.e ************ /n"));}
            else if (k == 18){mess.append(QString( "********** Betweenness Centrality Clustering ************ /n")); times = 1;}
            else if (k == 19){mess.append(QString( "********** Random Functional Digraph *************** /n"));}
            else{}
            for (int i = 0 ; i < times; i++)
            {
                if (k == 0) {continue;G.random_aggregate();}
                else if (k == 1) {continue;G.random_aggregate_with_degree_comparison(); }
                else if (k == 2) {continue;G.random_aggregate_with_weight_comparison(); }
                else if (k == 3) {continue;G.random_aggregate_with_neighbour_initial_degree_bias(); }
                else if (k == 4){continue;G.random_aggregate_with_neighbour_initial_degree_bias_with_constraint();}
                else if (k == 5) {continue;G.random_aggregate_with_neighbour_CURRENT_degree_bias(); }
                else if (k == 6){continue;G.random_aggregate_with_neighbour_CURRENT_degree_bias_with_constraint();}
                else if (k == 7) {continue;G.random_aggregate_highest_CURRENT_degree_neighbour();}
                else if (k == 8) {continue;G.random_aggregate_with_minimum_weight_neighbour();}
                else if (k == 9) {continue;G.random_aggregate_probabilistic_lowest_degree_neighbour_destructive();}
                else if (k == 10) {continue;G.random_aggregate_probabilistic_candidate_with_minimum_weight_neighbour();}
                else if (k == 11) {continue;G.random_aggregate_greedy_max_degree();}
                else if (k == 12){continue;G.random_aggregate_greedy_max_weight();}
                else if (k == 13){continue;G.random_aggregate_retain_vertex_using_triangulation();}
                else if (k == 14){continue;G.random_aggregate_retain_vertex_using_probabilistic_triangulation();}
                else if (k == 15){continue;G.random_aggregate_with_highest_triangulated_vertex();}
                else if (k == 16){continue;G.random_aggregate_retain_vertex_using_triangulation_times_weight();}
                else if (k == 17){continue;G.random_aggregate_retain_vertex_using_triangulation_of_cluster();}
                else if (k == 18){continue;G.betweenness_centrality_clustering();}
                else if (k == 19){G.random_functional_digraph();}
                double q = G.LARGE_compute_modularity();
                iQ += q;
                //count cc
                iCC += G.count_result_connected_component();
            }
            //add to list
            double  normalisedQ = iQ/times,
                    normalisedCC = (double) iCC/times;
            sQ << normalisedQ; sCC << normalisedCC;
            //write to file
            file.open(QFile::Text | QFile::Append);
            QTextStream out(&file);
            mess.append(QString("%1/t%2/t%3/t%4/t%5/t%6/n").arg(normalisedQ).arg(normalisedCC).arg(n));
            out << mess;
            file.close();
        }
        Q << sQ; CC << sCC;
    }
    writeSeperateFile(QString("Q"), Q);
  //  plotMatrix(Q, QString("Q"));
    qDebug() << CC;
    plotIntMatrix(CC, QString("Connected_Component"));
}



/////////////////////////////////////////////////////////////////////////////////////////
/// /brief Experiment with K Circulant Graph
/// /param k
///
void generate_K_circulant(const int &k)
{
    int n = 100;
    TRnd rand;
    PUNGraph TGraph = PUNGraph::New();
    TGraph = TSnap::GenSmallWorld(n, k, 0.3, rand);
    TSnap::SaveEdgeList(TGraph,
                        "C:/Users/Dumex/Desktop/SocialNetworksCollection/K_Circulant/edge_file.txt",
                        "K Circulant");
}

void K_circulant_experiment()
{
    Graph G;
    QFile file("C:/Users/Dumex/Desktop/SocialNetworksCollection/GirvanNewmanExperiment/Stat.txt");
    QList<QList<double> > Q;
    QList<QList<int> > CC;
    for (int NodeOutDeg = 5 ; NodeOutDeg <= 50; NodeOutDeg+=5)
    {
        QList<double> sQ; //s = Set
        QList<int> sCC;
        //generate Watts Schorgat
        generate_K_circulant(NodeOutDeg);
        G.LARGE_hard_reset();
        G.read_edge(QString("C:/Users/Dumex/Desktop/SocialNetworksCollection/K_Circulant/"));
        for (int k = 0; k <= 18; k++)
        {
            //generate GN graph
            //run algorithm here
            int times = 10;
            double iQ = 0.0; // i = iterator
            quint32 iCC = 0;
            QString mess;
            if (k == 0) {mess.append(QString("********** I.a ************ /n"));}
            else if (k == 1) {mess.append(QString( "********** I.b ************ /n"));}
            else if (k == 2) {mess.append(QString( "********** I.c ************ /n"));}
            else if (k == 3) {mess.append(QString( "********** II.a ************ /n"));}
            else if (k == 4){mess.append(QString( "********** II.a(i) ************ /n"));}
            else if (k == 5) {mess.append(QString( "********** II.b ************ /n"));}
            else if (k == 6){mess.append(QString( "********** II.b(i) ************ /n"));}
            else if (k == 7) {mess.append(QString( "********** II.c ************ /n"));}
            else if (k == 8) {mess.append(QString( "********** II.d ************ /n"));}
            else if (k == 9) {mess.append(QString( "********** II.e ************ /n"));}
            else if (k == 10) {mess.append(QString( "********** II.f ************ /n"));}
            else if (k == 11) {mess.append(QString( "********** II.g ************ /n"));}
            else if (k == 12){mess.append(QString( "********** II.h ************ /n"));}
            else if (k == 13){mess.append(QString( "********** III.a ************ /n"));}
            else if (k == 14){mess.append(QString( "********** III.b ************ /n"));}
            else if (k == 15){mess.append(QString( "********** III.c ************ /n"));}
            else if (k == 16){mess.append(QString( "********** III.d ************ /n"));}
            else if (k == 17){mess.append(QString( "********** III.e ************ /n"));}
            else if (k == 18){mess.append(QString( "********** Betweenness Centrality Clustering ************ /n")); times = 1;}
            else{}
            for (int i = 0 ; i < times; i++)
            {
                if (k == 0) {G.random_aggregate();}
                else if (k == 1) {G.random_aggregate_with_degree_comparison(); }
                else if (k == 2) {G.random_aggregate_with_weight_comparison(); }
                else if (k == 3) {G.random_aggregate_with_neighbour_initial_degree_bias(); }
                else if (k == 4){G.random_aggregate_with_neighbour_initial_degree_bias_with_constraint();}
                else if (k == 5) {G.random_aggregate_with_neighbour_CURRENT_degree_bias(); }
                else if (k == 6){G.random_aggregate_with_neighbour_CURRENT_degree_bias_with_constraint();}
                else if (k == 7) {G.random_aggregate_highest_CURRENT_degree_neighbour();}
                else if (k == 8) {G.random_aggregate_with_minimum_weight_neighbour();}
                else if (k == 9) {G.random_aggregate_probabilistic_lowest_degree_neighbour_destructive();}
                else if (k == 10) {G.random_aggregate_probabilistic_candidate_with_minimum_weight_neighbour();}
                else if (k == 11) {G.random_aggregate_greedy_max_degree();}
                else if (k == 12){G.random_aggregate_greedy_max_weight();}
                else if (k == 13){G.random_aggregate_retain_vertex_using_triangulation();}
                else if (k == 14){G.random_aggregate_retain_vertex_using_probabilistic_triangulation();}
                else if (k == 15){G.random_aggregate_with_highest_triangulated_vertex();}
                else if (k == 16){G.random_aggregate_retain_vertex_using_triangulation_times_weight();}
                else if (k == 17){G.random_aggregate_retain_vertex_using_triangulation_of_cluster();}
                else if (k == 18){continue;G.betweenness_centrality_clustering();}
                double q = G.LARGE_compute_modularity();
                iQ += q;
                //count cc
                iCC += G.count_result_connected_component();
            }
            //add to list
            double  normalisedQ = iQ/times,
                    normalisedCC = (double) iCC/times;
            sQ << normalisedQ; sCC << normalisedCC;
            //write to file
            file.open(QFile::Text | QFile::Append);
            QTextStream out(&file);
            mess.append(QString("%1/t%2/t%3/t%4/t%5/t%6/n").arg(normalisedQ).arg(normalisedCC).arg(NodeOutDeg));
            out << mess;
            file.close();
        }
        Q << sQ; CC << sCC;
    }
    writeSeperateFile(QString("Q"), Q);
    plotMatrix(Q, QString("Q"));
    qDebug() << CC;
    plotIntMatrix(CC, QString("Connected_Component"));
}


//////////////////////////////////////////////////////////////////////////////////////////
/// /brief Experiment with Simple Cycle
/// /param n
/// /return
///
int RandomMappingOnACycle(const int& n)
{
    std::random_device rd;
    std::mt19937 generator(rd());
    generator.seed(qSqrt(QTime::currentTime().msec()));
    std::uniform_real_distribution<double> dis;
    QList<QPair<int,int> > edge;
    for(int i = 0; i < n; i++)
    {
        double ran = dis(generator);
        int u = -1;
        if (i == 0)
        {
            if (ran <= 0.5)  u = i+1;
            else            u = n-1;
        }
        else if (i == n-1)
        {
            if (ran <= 0.5)  u = 0;
            else            u = n-1;
        }
        else
        {
            if (ran <= 0.5)  u = i+1;
            else            u = i-1;
        }
        if (u == -1){ qDebug() << "ERROR WHILE GENERATING EDGE"; return -1;}
        edge.append(qMakePair(i,u));
    }
    //check number of conencted component
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> NormalGraph;
    NormalGraph G;
    for (int i = 0; i < n; i++)
        boost::add_vertex(G);
    for (int i = 0; i < edge.size(); i++)
    {
        QPair<int,int> p = edge.at(i);
        boost::add_edge(p.first, p.second, G);
    }
    std::vector<int> component(boost::num_vertices(G));
    int num = boost::connected_components(G, &component[0]);
    return num;
}


void RandomMappingOnACycle_Exp()
{
    Graph G;
    G.manual_set_working_dir(workingDir);
    int n = 10;
    TVec<TPair<TInt, TInt> >    I_a, I_b, I_c,
            II_a, II_a_i, II_b, II_b_i, II_c, II_d, II_e, II_f, II_g, II_h,
            R1a,
            I_x,
            RFD;
    while (n <= 10000)
    {
        G.generateSimpleCycle(n);
        writeEdgeFile(G);
        int times = 10;
        for (int k = 0; k <= 15; k++)
        {
            int sum_cc = 0;
            for (int i = 0 ; i < times; i++)
            {
                int cc = 0;
                if (k == 0) {G.random_aggregate();}
                else if (k == 1) {G.random_aggregate_with_degree_comparison(); }
                else if (k == 2) {G.random_aggregate_with_weight_comparison(); }
                else if (k == 3) {G.random_aggregate_with_neighbour_initial_degree_bias(); }
                else if (k == 4){G.random_aggregate_with_neighbour_initial_degree_bias_with_constraint();}
                else if (k == 5) {G.random_aggregate_with_neighbour_CURRENT_degree_bias(); }
                else if (k == 6){G.random_aggregate_with_neighbour_CURRENT_degree_bias_with_constraint();}
                else if (k == 7) {G.random_aggregate_highest_CURRENT_degree_neighbour();}
                else if (k == 8) {G.random_aggregate_with_minimum_weight_neighbour();}
                else if (k == 9) {G.random_aggregate_probabilistic_lowest_degree_neighbour_destructive();}
                else if (k == 10) {G.random_aggregate_probabilistic_candidate_with_minimum_weight_neighbour();}
                else if (k == 11) {G.random_aggregate_greedy_max_weight();}
                else if (k == 12){G.random_aggregate_greedy_max_degree();}
                else if (k == 13){G.reverse_random_aggregate();}
                else if (k == 14){G.reverse_random_aggregate_with_degree_comparison();}
                else if (k == 15){cc = RandomMappingOnACycle(n); //RFD
                    if (cc == -1){qDebug() << "ERROR!, Terminating ..."; return;}
                    else{sum_cc += cc;continue;}
                }
                cc = G.count_result_connected_component();
                sum_cc += cc;
                G.LARGE_reset();
            }
            //add result here
            int normalised_cc = sum_cc / times;
            TPair<TInt, TInt> p(n,normalised_cc);
            if (k == 0) {I_a.Add(p);}
            else if (k == 1) {I_b.Add(p); }
            else if (k == 2) {I_c.Add(p);}
            else if (k == 3) {II_a.Add(p); }
            else if (k == 4){II_a_i.Add(p);}
            else if (k == 5) {II_b.Add(p); }
            else if (k == 6){II_b_i.Add(p);}
            else if (k == 7) {II_c.Add(p);}
            else if (k == 8) {II_d.Add(p);}
            else if (k == 9) {II_e.Add(p);}
            else if (k == 10) {II_f.Add(p);}
            else if (k == 11) {II_h.Add(p);}
            else if (k == 12){II_g.Add(p);}
            else if (k == 13){R1a.Add(p);}
            else if (k == 14){I_x.Add(p);}
            else if (k == 15){RFD.Add(p);} //RFD
        }
        G.LARGE_hard_reset();
        n+=50;
    }
    TGnuPlot Gp("RandomMappingOnACycleEXP", "Number of Connected Component As a Function of n");
    Gp.AddPlot(I_a,gpwLinesPoints, "I.a");
    Gp.AddPlot(I_b,gpwLinesPoints, "I.b");
    Gp.AddPlot(I_c,gpwLinesPoints, "I.c");
    Gp.AddPlot(II_a,gpwLinesPoints, "II.a");
    Gp.AddPlot(II_a_i,gpwLinesPoints, "II.a.i");
    Gp.AddPlot(II_b,gpwLinesPoints, "II.b");
    Gp.AddPlot(II_b_i,gpwLinesPoints, "II.b.i");
    Gp.AddPlot(II_c,gpwLinesPoints, "II.c");
    Gp.AddPlot(II_d,gpwLinesPoints, "II.d");
    Gp.AddPlot(II_e,gpwLinesPoints, "II.e");
    Gp.AddPlot(II_f,gpwLinesPoints, "II.f");
    Gp.AddPlot(II_g,gpwLinesPoints, "II.g");
    Gp.AddPlot(II_h,gpwLinesPoints, "II.h");
    Gp.AddPlot(R1a,gpwLinesPoints, "R1A");
    Gp.AddPlot(I_x,gpwLinesPoints, "I.x");
    Gp.AddPlot(RFD,gpwLinesPoints, "RFD");
    Gp.SetXYLabel("n", "Connected Component");
    Gp.SavePng();
}

///////////////////////////////////////////////////////////////////////////////////////////
/// /brief Experiment with Binary Tree
/// /param h
/// /return
///
int RandomMappingOnBinaryTree(const int &h)
{
    unsigned int n = qPow(2,h)-1;
    std::random_device rd;
    std::mt19937 generator(rd());
    generator.seed(qSqrt(QTime::currentTime().msec()));
    std::uniform_real_distribution<double> dis;
    QList<QPair<int,int> > edge;
    for(int i = 0; i < n; i++)
    {
        int level = qFloor(std::log2(i+1))+1;
        double ran = dis(generator);
        unsigned int parent = qFloor((i-1)/2),
                     leftC = i*2 + 1,
                     rightC = i*2 + 2,
                     chosen = -1;
        if (level == 1)
        {
            if (ran <= 0.5) chosen = leftC;
            else            chosen = rightC;
        }
        else if (level > 1 && level <= h-1)
        {
            if (ran <= 0.33333) chosen = leftC;
            else if (ran > 0.33333 && ran <= 0.66666)   chosen = rightC;
            else    chosen = parent;
        }
        else if (level == h)
        {
            chosen = parent;
        }
        if (chosen == -1){qDebug() << "Error While Generating Edge..."; return -1;}
        edge.append(qMakePair(i,chosen));

    }

    //check number of conencted component
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> NormalGraph;
    NormalGraph G;
    for (int i = 0; i < n; i++)
        boost::add_vertex(G);
    for (int i = 0; i < edge.size(); i++)
    {
        QPair<int,int> p = edge.at(i);
        boost::add_edge(p.first, p.second, G);
    }
    std::vector<int> component(boost::num_vertices(G));
    int num = boost::connected_components(G, &component[0]);
    qDebug() << edge;
    qDebug() << num;
    return -1;
    return num;
}


void RandomMappingOnBinaryTree_EXP()
{
    Graph G;
    G.manual_set_working_dir(workingDir);
    int h = 4;
    TVec<TPair<TInt, TInt> > RM, Ia, Ib, Ic, IIa, IIb,  IIc, IId, IIe, IIf, IIg, IIh;
    while (h <= 20)
    {
        int n = qPow(2,h)-1;
        int cc = RandomMappingOnBinaryTree(h); //random mapping
        if (cc == -1) return;
        TPair<TInt,TInt> p(n,cc);
        RM.Add(p);
        G.generateBinaryTree(h);
        writeEdgeFile(G);
        for(int i = 0; i <= 10; i++)
        {
            if (i == 0)
            {
                G.random_aggregate();
                int c = G.count_result_connected_component();
                TPair<TInt,TInt> p2(n,c);
                Ia.Add(p2);
            }
            else if (i == 1)
            {
                G.random_aggregate_with_degree_comparison();
                int c = G.count_result_connected_component();
                TPair<TInt,TInt> p2(n,c);
                Ib.Add(p2);
            }
            else if (i == 2)
            {
                G.random_aggregate_with_weight_comparison();
                int c = G.count_result_connected_component();
                TPair<TInt,TInt> p2(n,c);
                Ic.Add(p2);
            }
            else if (i == 3) {
                G.random_aggregate_with_neighbour_initial_degree_bias();
                int c = G.count_result_connected_component();
                TPair<TInt,TInt> p2(n,c);
                IIa.Add(p2);}
            else if (i == 4) {G.random_aggregate_with_neighbour_CURRENT_degree_bias();
                int c = G.count_result_connected_component();
                TPair<TInt,TInt> p2(n,c);
                IIb.Add(p2);}
            else if (i == 5) {G.random_aggregate_highest_CURRENT_degree_neighbour();
                int c = G.count_result_connected_component();
                TPair<TInt,TInt> p2(n,c);
                IIc.Add(p2);}
            else if (i == 6) {G.random_aggregate_with_minimum_weight_neighbour();
                int c = G.count_result_connected_component();
                TPair<TInt,TInt> p2(n,c);
                IId.Add(p2);}
            else if (i == 7) {G.random_aggregate_probabilistic_lowest_degree_neighbour_destructive();
                int c = G.count_result_connected_component();
                TPair<TInt,TInt> p2(n,c);
                IIe.Add(p2);}
            else if (i == 8) {G.random_aggregate_probabilistic_candidate_with_minimum_weight_neighbour();
                int c = G.count_result_connected_component();
                TPair<TInt,TInt> p2(n,c);
                IIf.Add(p2);}
            else if (i == 9) {G.random_aggregate_greedy_max_weight();
                int c = G.count_result_connected_component();
                TPair<TInt,TInt> p2(n,c);
                IIg.Add(p2);}
            else if (i == 10){G.random_aggregate_greedy_max_degree();
                int c = G.count_result_connected_component();
                TPair<TInt,TInt> p2(n,c);
                IIh.Add(p2);}
        }
        G.LARGE_hard_reset();
        h++;
    }
    TGnuPlot Gp("RandomMappingOnBinaryTreeEXP", "Number of Connected Component As a Function of n");
    Gp.AddPlot(RM,gpwLinesPoints, "Random Mapping");
    Gp.AddPlot(Ia,gpwLinesPoints, "I.a");
    Gp.AddPlot(Ib,gpwLinesPoints, "I.b");
    Gp.AddPlot(Ic,gpwLinesPoints, "I.c");
    Gp.AddPlot(IIa,gpwLinesPoints, "II.a");
    Gp.AddPlot(IIb,gpwLinesPoints, "II.b");
    Gp.AddPlot(IIc,gpwLinesPoints, "II.c");
    Gp.AddPlot(IId,gpwLinesPoints, "II.d");
    Gp.AddPlot(IIe,gpwLinesPoints, "II.e");
    Gp.AddPlot(IIf,gpwLinesPoints, "II.f");
    Gp.AddPlot(IIg,gpwLinesPoints, "II.g");
    Gp.AddPlot(IIh,gpwLinesPoints, "II.h");
    Gp.SetXYLabel("n", "Connected Component");
    Gp.SavePng();
}

//calculate the coefficient given by the generating function of E(x)
// E(x) = A(x)B(x) = e^(x^2/2).e^(-x)
unsigned int factorial(int n)
{
    if (n == 0) return 1;
    unsigned int fac = 1;
    for(int i = 1; i <= n; i++)
    {
        fac *= i;
    }
    return fac;
}
double xnGx(int n)
{
    double coeff = 0.0;
    int i = 0;
    if (n%2 != 0)   i = 1;
    for(i; i <= n; i+=2)
    {
        int j = (n - i)/2;
        quint32 fac_i = factorial(i),
                fac_j = factorial(j),
                pow = qPow(2,j);
        double p = (double) 1.0/pow;
        coeff += (double) 1/fac_i*1/fac_j*p;
    }
    if (n%2 == 0) return coeff;
    else    return coeff*-1;
}
double xnEx(int n)
{
    double coeff = 0.0;
    for(int i = 0; i <= n; i++)
    {
       coeff += xnGx(i)*(n+1-i);
    }
    return coeff;
}



///////////////////////////////////////////////////////////////////////////////////////////
/// /brief RUN ABOVE EXPERIMENTS WITH PREDETERMINED GRAPHS
/// /param filepath
/// /return
///
void SocialGraph_exp(const QString &graph_name)
{
    QList<QList<double> > RAND, JACCARD, ARI, Q, GN;
    for (int t = 0 ; t < 1; t++)
    {
        QList<double> sRAND, sJACCARD, sARI, sQ, sGN; //s = Set
        for (int k = 0; k <= 22; k++)
        {
            //generate GN graph
            int times = 10;
            double iRAND = 0.0 , iJACCARD = 0.0, iARI = 0.0, iQ = 0.0, iGN = 0.0; // i = iterator
            //run algorithm here
            QString mess;
            if (k == 0) {mess.append(QString("********** I.a ************ /n"));}
            else if (k == 1) {mess.append(QString( "********** I.b ************ /n"));}
            else if (k == 2) {mess.append(QString( "********** I.c ************ /n"));}
            else if (k == 3) {mess.append(QString( "********** II.a ************ /n"));}
            else if (k == 4){mess.append(QString( "********** II.a(i) ************ /n"));}
            else if (k == 5) {mess.append(QString( "********** II.b ************ /n"));}
            else if (k == 6){mess.append(QString( "********** II.b(i) ************ /n"));}
            else if (k == 7) {mess.append(QString( "********** II.c ************ /n"));}
            else if (k == 8) {mess.append(QString( "********** II.d ************ /n"));}
            else if (k == 9) {mess.append(QString( "********** II.e ************ /n"));}
            else if (k == 10) {mess.append(QString( "********** II.f ************ /n"));}
            else if (k == 11) {mess.append(QString( "********** II.h ************ /n"));}
            else if (k == 12){mess.append(QString( "********** II.g ************ /n"));}
            else if (k == 13){mess.append(QString( "********** III.a ************ /n"));}
            else if (k == 14){mess.append(QString( "********** III.b ************ /n"));}
            else if (k == 15){mess.append(QString( "********** III.c ************ /n"));}
            else if (k == 16){mess.append(QString( "********** III.d ************ /n"));}
            else if (k == 17){mess.append(QString( "********** III.e ************ /n"));}
            else if (k == 18){mess.append(QString( "********** Betweenness Centrality Clustering ************ /n")); times = 1;}
            else if (k == 19){mess.append(QString( "********** CNM Clustering ************ /n")); times = 1;}
            else if (k == 20){mess.append(QString( "********** R1a ************ /n"));}
            else if (k == 21){mess.append(QString( "********** I.x ************ /n"));}
            else if (k == 22){mess.append(QString( "********** RFD ************ /n"));}
            else{}
            for (int i = 0 ; i < times; i++)
            {
                Graph G;
                G.manual_set_working_dir(workingDir);
                G.read_edge(workingDir);
                G.load_ground_truth_communities();
                if (k == 0) {G.random_aggregate();}
                else if (k == 1) {G.random_aggregate_with_degree_comparison(); }
                else if (k == 2) {G.random_aggregate_with_weight_comparison(); }
                else if (k == 3) {G.random_aggregate_with_neighbour_initial_degree_bias(); }
                else if (k == 4){G.random_aggregate_with_neighbour_initial_degree_bias_with_constraint();}
                else if (k == 5) {G.random_aggregate_with_neighbour_CURRENT_degree_bias(); }
                else if (k == 6){G.random_aggregate_with_neighbour_CURRENT_degree_bias_with_constraint();}
                else if (k == 7) {G.random_aggregate_highest_CURRENT_degree_neighbour();}
                else if (k == 8) {G.random_aggregate_with_minimum_weight_neighbour();}
                else if (k == 9) {G.random_aggregate_probabilistic_lowest_degree_neighbour_destructive();}
                else if (k == 10) {G.random_aggregate_probabilistic_candidate_with_minimum_weight_neighbour();}
                else if (k == 11) {G.random_aggregate_greedy_max_weight();}
                else if (k == 12){G.random_aggregate_greedy_max_degree();}
                else if (k == 13){G.random_aggregate_retain_vertex_using_triangulation();}
                else if (k == 14){G.random_aggregate_retain_vertex_using_probabilistic_triangulation();}
                else if (k == 15){G.random_aggregate_with_highest_triangulated_vertex();}
                else if (k == 16){G.random_aggregate_retain_vertex_using_triangulation_times_weight();}
                else if (k == 17){G.random_aggregate_retain_vertex_using_triangulation_of_cluster();}
                else if (k == 18){G.betweenness_centrality_clustering();}
                else if (k == 19){G.fast_CMN();}
                else if (k == 20){G.reverse_random_aggregate();}
                else if (k == 21){G.reverse_random_aggregate_with_degree_comparison();}
                else if (k == 22){G.random_functional_digraph();}
                QList<double> id = G.LARGE_compute_Pairwise_efficient(-1);
                iRAND+=id[0];
                iJACCARD+=id[1];
                iARI+=id[2];
                iQ+=G.LARGE_compute_modularity();
               // double gn = G.compute_GN_index();
                double gn = 0.0;
                G.LARGE_hard_reset();
                iGN += gn;
            }
            //add to list
            double normalisedRand = iRAND/times,
                    normalisedJaccard = iJACCARD/times,
                    normalisedARI = iARI/times,
                    normalisedQ = iQ/times,
                    normalisedGN = (double) iGN/times;
            sRAND << normalisedRand; sJACCARD << normalisedJaccard; sARI << normalisedARI; sQ << normalisedQ; sGN << normalisedGN;
            //write to file
            QFile file(workingDir +"Stat.txt");
            file.open(QFile::Text | QFile::Append);
            QTextStream out(&file);
            mess.append(QString("%1/t%2/t%3/t%4/t%5/t%6/n").arg(normalisedRand).arg(normalisedJaccard).arg(normalisedARI).arg(normalisedQ).arg(normalisedGN).arg(t));
            out << mess;
            file.close();
        }
        RAND << sRAND; JACCARD << sJACCARD; ARI << sARI; Q << sQ; GN << sGN;
    }
    writeSeperateFile(QString("ARI"), ARI);
    writeSeperateFile(QString("JACCARD"), JACCARD);
    writeSeperateFile(QString("Q"), Q);
    writeSeperateFile(QString("GN"), GN);
    plotMatrix(ARI, QString("ARI"));
    plotMatrix(JACCARD, QString("Jaccard"));
    plotMatrix(Q, QString("Q"));
    plotMatrix(GN, QString("FoCC"));
}


///////////////////////////////////////////////////////////////////////////////////////////
/// /brief RUN ABOVE EXPERIMENTS WITH PREDETERMINED GRAPHS
/// /param filepath
/// /return
///
void LFR_benchmark_exp()
{
    qDebug() << "LFR Experiment Started";
    QList<int> blocked_type;
    blocked_type << 9 << 10 << 14 << 15 << 16 << 18 << 19//eliminated after n > 10,000
                 << 0 << 1 << 2 << 3 << 4 << 5 << 6 << 7 << 8 << 19 << 11 << 12; //eliminated after n > 100,000

    QList<QList<double> > RAND, JACCARD, ARI, Q, GN;
    for (int t = 0 ; t < 1; t++)
    {
        QList<double> sRAND, sJACCARD, sARI, sQ, sGN; //s = Set
        for (int k = 0; k <= 22; k++)
        {
            //generate GN graph
            int times = 5;
            double iRAND = 0.0 , iJACCARD = 0.0, iARI = 0.0, iQ = 0.0, iGN = 0.0; // i = iterator
            //run algorithm here
            QString mess;
            if (k == 0) {mess.append(QString("********** I.a ************ /n"));}
            else if (k == 1) {mess.append(QString( "********** I.b ************ /n"));}
            else if (k == 2) {mess.append(QString( "********** I.c ************ /n"));}
            else if (k == 3) {mess.append(QString( "********** II.a ************ /n"));}
            else if (k == 4){mess.append(QString( "********** II.a(i) ************ /n"));}
            else if (k == 5) {mess.append(QString( "********** II.b ************ /n"));}
            else if (k == 6){mess.append(QString( "********** II.b(i) ************ /n"));}
            else if (k == 7) {mess.append(QString( "********** II.c ************ /n"));}
            else if (k == 8) {mess.append(QString( "********** II.d ************ /n"));}
            else if (k == 9) {mess.append(QString( "********** II.e ************ /n"));}
            else if (k == 10) {mess.append(QString( "********** II.f ************ /n"));}
            else if (k == 11) {mess.append(QString( "********** II.h ************ /n"));}
            else if (k == 12){mess.append(QString( "********** II.g ************ /n"));}
            else if (k == 13){mess.append(QString( "********** III.a ************ /n"));}
            else if (k == 14){mess.append(QString( "********** III.b ************ /n"));}
            else if (k == 15){mess.append(QString( "********** III.c ************ /n"));}
            else if (k == 16){mess.append(QString( "********** III.d ************ /n"));}
            else if (k == 17){mess.append(QString( "********** III.e ************ /n"));}
            else if (k == 18){mess.append(QString( "********** Betweenness Centrality Clustering ************ /n")); times = 1;}
            else if (k == 19){mess.append(QString( "********** CNM Clustering ************ /n")); times = 1;}
            else if (k == 20){mess.append(QString( "********** R1a ************ /n"));}
            else if (k == 21){mess.append(QString( "********** I.x ************ /n"));}
            else if (k == 22){mess.append(QString( "********** RFD ************ /n"));}
            else{}
            for (int i = 0 ; i < times; i++)
            {
                if (blocked_type.contains(k))
                {
                    iRAND+=0;
                    iJACCARD+=0;
                    iARI+=0;
                    iQ+=0;
                    iGN+=0;
                }
                else
                {
                    Graph G;
                    G.manual_set_working_dir(workingDir);
                    G.read_edge(workingDir);
                    G.load_LFR_groundTruth();
                    if (k == 0) {G.random_aggregate();}
                    else if (k == 1) {G.random_aggregate_with_degree_comparison(); }
                    else if (k == 2) {G.random_aggregate_with_weight_comparison(); }
                    else if (k == 3) {G.random_aggregate_with_neighbour_initial_degree_bias(); }
                    else if (k == 4){G.random_aggregate_with_neighbour_initial_degree_bias_with_constraint();}
                    else if (k == 5) {G.random_aggregate_with_neighbour_CURRENT_degree_bias(); }
                    else if (k == 6){G.random_aggregate_with_neighbour_CURRENT_degree_bias_with_constraint();}
                    else if (k == 7) {G.random_aggregate_highest_CURRENT_degree_neighbour();}
                    else if (k == 8) {G.random_aggregate_with_minimum_weight_neighbour();}
                    else if (k == 9) {G.random_aggregate_probabilistic_lowest_degree_neighbour_destructive();}
                    else if (k == 10) {G.random_aggregate_probabilistic_candidate_with_minimum_weight_neighbour();}
                    else if (k == 11) {G.random_aggregate_greedy_max_weight();}
                    else if (k == 12){G.random_aggregate_greedy_max_degree();}
                    else if (k == 13){G.random_aggregate_retain_vertex_using_triangulation();}
                    else if (k == 14){G.random_aggregate_retain_vertex_using_probabilistic_triangulation();}
                    else if (k == 15){G.random_aggregate_with_highest_triangulated_vertex();}
                    else if (k == 16){G.random_aggregate_retain_vertex_using_triangulation_times_weight();}
                    else if (k == 17){G.random_aggregate_retain_vertex_using_triangulation_of_cluster();}
                    else if (k == 18){G.betweenness_centrality_clustering();}
                    else if (k == 19){G.fast_CMN();}
                    else if (k == 20){G.reverse_random_aggregate();}
                    else if (k == 21){G.reverse_random_aggregate_with_degree_comparison();}
                    else if (k == 22){G.random_functional_digraph();}
                    QList<double> id = G.LARGE_compute_Pairwise_efficient(-1);
                    if (id[0] > 1 || id[1] > 1 || id[2] > 1)
                    {
                        qDebug() << "INCORRECT VALUES OF PAIRWISE COMPARISON";
                        return;
                    }
                    iRAND+=id[0];
                    iJACCARD+=id[1];
                    iARI+=id[2];
                    iQ+=G.LARGE_compute_modularity();
                   // double gn = G.compute_GN_index();
                    double gn = 0.0;
                    G.LARGE_hard_reset();
                    iGN += gn;
                }
            }
            //add to list
            double normalisedRand = iRAND/times,
                    normalisedJaccard = iJACCARD/times,
                    normalisedARI = iARI/times,
                    normalisedQ = iQ/times,
                    normalisedGN = (double) iGN/times;
            sRAND << normalisedRand; sJACCARD << normalisedJaccard; sARI << normalisedARI; sQ << normalisedQ; sGN << normalisedGN;
            //write to file
            QFile file(workingDir +"Stat.txt");
            file.open(QFile::Text | QFile::Append);
            QTextStream out(&file);
            mess.append(QString("%1 /t %2 /t %3 /t %4 /t %5 /t %6/n").arg(normalisedRand).arg(normalisedJaccard).arg(normalisedARI).arg(normalisedQ).arg(normalisedGN).arg(t));
            out << mess;
            file.close();
        }
        RAND << sRAND; JACCARD << sJACCARD; ARI << sARI; Q << sQ; GN << sGN;
    }
    writeSeperateFile(QString("ARI"), ARI);
   // writeSeperateFile(QString("JACCARD"), JACCARD);
    writeSeperateFile(QString("Q"), Q);
  //  writeSeperateFile(QString("GN"), GN);
    plotMatrix(ARI, QString(workingDir + "ARI"));
    plotMatrix(JACCARD, QString(workingDir + "Jaccard"));
    plotMatrix(Q, QString(workingDir + "Q"));
  //  plotMatrix(GN, QString("FoCC"));
}


int main(int argc, char *argv[])
{
   // qInstallMessageHandler(myMessageOutput);
    //Large Test for LFR
    qDebug() << "Large Gnp Experiments Initialised";
    workingDir = "C:/Users/Dumex/Desktop/SocialNetworksCollection/Triangle_Gnp/";
    LARGE_Gnp_experiment();
    return 0;

/*
    while (true)
    {
        std::cout << "Enter: ";
        int k;
        std::cin >> k;
        if (k == 0) {G.random_aggregate();}
        else if (k == 1) {G.random_aggregate_with_degree_comparison(); }
        else if (k == 2) {G.random_aggregate_with_weight_comparison(); }
        else if (k == 3) {G.random_aggregate_with_neighbour_initial_degree_bias(); }
        else if (k == 4){G.random_aggregate_with_neighbour_initial_degree_bias_with_constraint();}
        else if (k == 5) {G.random_aggregate_with_neighbour_CURRENT_degree_bias(); }
        else if (k == 6){G.random_aggregate_with_neighbour_CURRENT_degree_bias_with_constraint();}
        else if (k == 7) {G.random_aggregate_highest_CURRENT_degree_neighbour();}
        else if (k == 8) {G.random_aggregate_with_minimum_weight_neighbour();}
        else if (k == 9) {G.random_aggregate_probabilistic_lowest_degree_neighbour_destructive();}
        else if (k == 10) {G.random_aggregate_probabilistic_candidate_with_minimum_weight_neighbour();}
        else if (k == 11) {G.random_aggregate_greedy_max_weight();}
        else if (k == 12){G.random_aggregate_greedy_max_degree();}
        else if (k == 13){G.random_aggregate_retain_vertex_using_triangulation();}
        else if (k == 14){G.random_aggregate_retain_vertex_using_probabilistic_triangulation();}
        else if (k == 15){G.random_aggregate_with_highest_triangulated_vertex();}
        else if (k == 16){G.random_aggregate_retain_vertex_using_triangulation_times_weight();}
        else if (k == 17){G.random_aggregate_retain_vertex_using_triangulation_of_cluster();}
        else if (k == 18){G.betweenness_centrality_clustering();}
        else if (k == 19) break;
    }
    std::cout << "Terminating";*/
   //
    return 0;
}
