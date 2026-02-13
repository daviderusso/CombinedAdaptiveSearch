#include <stdio.h>
#include <stdlib.h>
#include "BinaryCombinedHeuristic.h"
#include <string.h>
#include <sys/stat.h>
#include <errno.h>


int main() {
    const char *filename = "data/W.xy"; // File name of the graph

    clock_t startReadGraph = clock();
    Graph **Grafi = readGraphAndCreateReverse(filename);

    Graph *graph = Grafi[0];
    Graph *reverseGraph = Grafi[1];
    double TimeReadGraph = (double) (clock() - startReadGraph) / CLOCKS_PER_SEC;

    // const char *filein = "data/160WInstancesHeur_AP.csv";
    const char *filein = "data/8WInstancesHeur_AP.csv";
    // const char *filein = "data/1WInstancesHeur_AP.csv";
    int TimeLimit = 60;

    const char *res_dir = "res";
    if (mkdir(res_dir, 0777) != 0) {
        if (errno != EEXIST) {
            perror("Errore creazione cartella res");
            return EXIT_FAILURE;
        }
    }

    char outfile[256];
    char filetime[256];
    snprintf(outfile, sizeof(outfile), "%s/results_heur_W.csv", res_dir);
    snprintf(filetime, sizeof(filetime), "%s/results_times_heur_W.csv", res_dir);
    FILE *fp = fopen(filein, "r");
    FILE *fout = fopen(outfile, "w");
    FILE *fouttime = fopen(filetime, "w");
    if (!fp || !fout) {
        perror("Errore apertura file");
        return EXIT_FAILURE;
    }

    char buffer[1024];
    // Skip the header
    if (!fgets(buffer, sizeof(buffer), fp)) {
        fprintf(stderr, "File vuoto o errore in lettura header\n");
        fclose(fp);
        return EXIT_FAILURE;
    }

    // Headers of the CSV
    fprintf(
        fout,
        "Source;Target;Budget;"
        "runtime_readGraph;runtime_SP_C1;runtime_SP_C2;"
        "SP1_C1;SP1_C2;SP2_C1;SP2_C2;"
        "runtime_Red;runtime_BCH;id_time_BCH;"
        "C1_BCH;C2_BCH;"
        "BestIter;BestLambda;NumIter\n");
    fprintf(fouttime,
            "Source;Target;Budget;"
            "runtime_readGraph;runtime_SPTreeSource;runtime_SPTreeDest;runtime_Red;"
            "runtime_first_lambda;runtime_BCH_total;"
            "perc_eliminated_nodes;perc_eliminated_arcs\n");

    char *line = buffer;

    int source, target;
    int W;
    // Read one line at a time
    while (fgets(buffer, sizeof(buffer), fp)) {
        line = buffer;
        // Removes '\n' in queue
        buffer[strcspn(buffer, "\r\n")] = 0;

        // Read data from the line: "%d;%d;%d"
        if (sscanf(line, "%d;%d;%d", &source, &target, &W) == 3) {
            printf("\nsource=%d  target=%d  Budget=%d\n", source, target, W);
        } else {
            fprintf(stderr, "Riga malformata: %s\n", line);
        }

        Results *Risultati = (Results *) malloc(sizeof(Results));

        clock_t startHeur = clock();
        source = source - 1; //per allineare gli id ai test lanciati nell'esatto
        target = target - 1; //per allineare gli id ai test lanciati nell'esatto
        BinaryCombinedHeuristic_AStar(graph, reverseGraph, source--, target--, W, TimeLimit, Risultati);
        double TimeTotalHeur = (double) (clock() - startHeur) / CLOCKS_PER_SEC;

        printf("\nIl cammino minimo consuma %lf risorse, e la sua lunghezza è %lf.", Risultati->SP1_C2,
               Risultati->SP1_C1);
        printf(
            "\nIl cammino Lambda consuma %lf risorse, e la sua lunghezza è %lf. Tempo individuazione: %lf s. Tempo impiegato: %lf s.\n",
            Risultati->BCH_C2, Risultati->BCH_C1, Risultati->idTime_BCH, Risultati->runtime_BCH);

        fprintf(fout, "%d;%d;%d;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf,%d;%f;%d\n",
                source, target, W,
                TimeReadGraph, Risultati->runtime_SP1, Risultati->runtime_SP2,
                Risultati->SP1_C1, Risultati->SP1_C2, Risultati->SP2_C1, Risultati->SP2_C2,
                Risultati->runtime_red, Risultati->runtime_BCH, Risultati->idTime_BCH,
                Risultati->BCH_C1, Risultati->BCH_C2,
                Risultati->bestIter, Risultati->bestLambda, Risultati->numIter);
        fflush(fout);

        fprintf(fouttime, "%d;%d;%d;%lf;%lf;%lf;%lf;%lf;%lf;%f;%f;\n",
                source, target, W,
                TimeReadGraph, Risultati->runtime_sptree_source, Risultati->runtime_sptree_destination,
                Risultati->runtime_red,
                Risultati->runtime_ALambda, TimeTotalHeur,
                Risultati->PercRedNodes, Risultati->PercRedArcs);
        fflush(fouttime);
        {
            char detailsFile[256];
            snprintf(detailsFile, sizeof(detailsFile), "%s/%d-%d-%d_sol_details.txt", res_dir, source, target, W);
            FILE *fdetails = fopen(detailsFile, "w");
            if (fdetails != NULL) {
                fprintf(fdetails, "c1 c2 timeFound lambda iter\n");
                for (int i = 0; i < Risultati->numImprovements; i++) {
                    const Improvement *imp = &Risultati->improvements[i];
                    fprintf(
                        fdetails,
                        "%d %d %lf %f %d\n",
                        (int) imp->C1Length,
                        (int) imp->C2Length,
                        imp->timeFound,
                        imp->lambda,
                        imp->iter
                    );
                }
                fclose(fdetails);
            } else {
                perror("Errore apertura file dettagli soluzioni");
            }
        }
        free(Risultati->improvements);
        free(Risultati);
        printf("\n---------------------------------\n");
    }

    fclose(fp);
    fclose(fout);
    freeGraph(graph);
    freeGraph(reverseGraph);
    return 0;
}

// IL GRAFO HA GLI ID CHE PARTONO DA 0 sia per noi che per l'esatto.
// per le istanze va fatto -1 ai nodi source e target
