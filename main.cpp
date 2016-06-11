#include <iostream>

#include <QCoreApplication>
#include <QtGlobal>
#include <QtDebug>
#include <QFile>
#include <QTextStream>

#include <stdio.h>
#include <stdlib.h>
#include "mygraph.h"

#include "Snap.h"


enum RandomAgg { I_a, I_b, I_c,
                 II_a, II_a_i, II_b, II_b_i, II_c, II_d, II_e, II_f, II_g, II_h,
                 III_a, III_b, III_c, III_d, III_e,
                 GN_Clustering
               };
const char *name[] = { "I.a", "I.b", "I_c",
                       "II.a", "II.a.i", "II.b", "II.b.i", "II.c", "II.d", "II.e", "II.f", "II.g", "II.h",
                       "III.a", "III.b", "III.c"," III.d", "III.e",
                       "GN_Clustering"
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
    TSnap::SaveEdgeList(TGraph,
                        "C:/Users/Dumex/Desktop/SocialNetworksCollection/GirvanNewmanExperiment/edge_file.txt",
                        "Save as tab-separated list of edges");
}

void writeSeperateFile(QString fileName, const QList<QList<double> > &matrix)
{
    qDebug() << "Writing File";
    QFile file("C:/Users/Dumex/Desktop/SocialNetworksCollection/GirvanNewmanExperiment/" + fileName + ".txt");
    file.open(QFile::WriteOnly | QFile::Text);
    QTextStream out(&file);
    out.setRealNumberPrecision(4);
    int h = matrix.size(), w = matrix[0].size();
    for (int i = 0; i < h; i++)
    {
        QString str;
        for (int j = 0; j < w; j++)
        {
            str.append(QString::number(matrix[i][j])).append('\t');
        }
        out << str;
        out << '\n';
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
            TPair<TInt, TFlt> p(i, matrix[i][j]);
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
        plotV.Add(Gp);
    }
    //plot
    for(int i = 0; i < AllSeries.Len(); i++)
    {
        if (i <= RandomAgg::I_c) plotV[0].AddPlot(AllSeries[i], gpwLinesPoints, name[i]);
        else if (i > RandomAgg::I_c && i <= RandomAgg::II_h) plotV[1].AddPlot(AllSeries[i], gpwLinesPoints, name[i]);
        else if (i > RandomAgg::II_h && i <= RandomAgg::III_e) plotV[2].AddPlot(AllSeries[i], gpwLinesPoints, name[i]);
        else
        {
            plotV[3].AddPlot(AllSeries[RandomAgg::I_c], gpwLinesPoints, name[RandomAgg::I_c]);
            plotV[3].AddPlot(AllSeries[RandomAgg::II_g], gpwLinesPoints, name[RandomAgg::II_g]);
            plotV[3].AddPlot(AllSeries[RandomAgg::III_a], gpwLinesPoints, name[RandomAgg::III_a]);
            plotV[3].AddPlot(AllSeries[RandomAgg::III_e], gpwLinesPoints, name[RandomAgg::III_e]);
            plotV[3].AddPlot(AllSeries[RandomAgg::GN_Clustering], gpwLinesPoints, name[RandomAgg::GN_Clustering]);
        }
    }
    //plot the best

    //export as png
    for (int i = 0; i < plotV.Len();i++)
    {
        TGnuPlot Gp = plotV[i];
        Gp.SetXYLabel("Z_out", label);
        Gp.SavePng();
    }
}

/** Girvan and Newman Experiments on 4 hidden partition of Gnp
 * @brief GN_experiment
 */
void GN_experiment()
{
    Graph G;
    QFile file("C:/Users/Dumex/Desktop/SocialNetworksCollection/GirvanNewmanExperiment/Stat.txt");
    QList<QList<double> > RAND, JACCARD, ARI, Q;
    for (int z_out = 0 ; z_out <= 12; z_out++)
    {
        QList<double> sRAND, sJACCARD, sARI, sQ; //s = Set
        for (int k = 0; k <= 18; k++)
        {
            //generate GN graph
            //run algorithm here
            double p_in = 0.0, p_out = 0.0;
            p_out = (double) z_out/96;
            p_in = (double) (16-z_out)/31;
            int times = 100;
            double iRAND = 0.0 , iJACCARD = 0.0, iARI = 0.0, iQ = 0.0; // i = iterator
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
            else{}
            for (int i = 0 ; i < times; i++)
            {
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
                QList<double> id = G.LARGE_compute_Pairwise_efficient(-1);
                iRAND+=id[0];
                iJACCARD+=id[1];
                iARI+=id[2];
                iQ+=G.LARGE_compute_modularity();
                G.LARGE_reset();
            }
            //add to list
            double normalisedRand = iRAND/times,
                    normalisedJaccard = iJACCARD/times,
                    normalisedARI = iARI/times,
                    normalisedQ = iQ/times;
            sRAND << normalisedRand; sJACCARD << normalisedJaccard; sARI << normalisedARI; sQ << normalisedQ;
            //write to file
            file.open(QFile::Text | QFile::Append);
            QTextStream out(&file);
            mess.append(QString("%1\t%2\t%3\t%4\t%5\n").arg(normalisedRand).arg(normalisedJaccard).arg(normalisedARI).arg(normalisedQ).arg(z_out));
            out << mess;
            file.close();
        }
        RAND << sRAND; JACCARD << sJACCARD; ARI << sARI; Q << sQ;
    }
    writeSeperateFile(QString("ARI"), ARI);
    writeSeperateFile(QString("JACCARD"), JACCARD);
    writeSeperateFile(QString("Q"), Q);
    plotMatrix(ARI, QString("ARI"));
    plotMatrix(JACCARD, QString("Jaccard"));
    plotMatrix(Q, QString("Q"));
}




int main(int argc, char *argv[])
{
   // qInstallMessageHandler(myMessageOutput);
    QString dirPath("C:/Users/Dumex/Desktop/SocialNetworksCollection/GirvanNewmanExperiment/");
    Graph G;
    G.manual_set_working_dir(dirPath);
    GN_experiment();

    qDebug() << "- Terminating ...";
    return 0;
}
