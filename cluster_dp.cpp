#include <cstdio>
#include <algorithm>
#include <cmath>

#include <iostream>

//#define IS_DEBUG
#ifdef IS_DEBUG
    #define DEBUG(x) cout << (#x) << " = " << (x) << endl
    #define RUN(x) x
#else
    #define DEBUG(x)
    #define RUN(x)
#endif

using namespace std;

const int MAX_LINES = 2000 * 2000;

struct SLink
{
    int u, v;
    double d;
};

SLink link[MAX_LINES];
int nlines;
int n;
double **dist;

double dc;
double *rho;
double *delta;
int *neigh;
int *rhoIdx;

int *cl;
int *icl;
int *halo;
int ncluster;

inline bool cmpLink(const SLink &a, const SLink &b)
{
    return a.d < b.d;
}

inline bool cmpRhoIdx(const int &a, const int &b)
{
    return rho[a] < rho[b];
}

void init()
{
    char filename[50];
    printf("The only input needed is a distance matrix file\n");
    printf("The format of this file should be: \n");
    printf("Column 1: id of element i\n");
    printf("Column 2: id of element j\n");
    printf("Column 3: dist(i,j)\n");
    printf("name of the distance matrix file?\n");
    scanf("%s", filename);
    RUN(strcpy(filename, "example_distances.dat"));
    
    printf("Reading input distance matrix\n");
    FILE *fp = fopen(filename, "r");
    nlines = 0;
    while (fscanf(fp, "%d %d %lf", &link[nlines].u, &link[nlines].v, &link[nlines].d) != EOF) {
        --link[nlines].u;
        --link[nlines].v;
        ++nlines;
    }
    fclose(fp);
    DEBUG(nlines);
    
    n = 0;
    for (int i = 0; i < nlines; ++i) {
        n = max(n, max(link[i].u, link[i].v));
    }
    ++n;
    DEBUG(n);
    
    dist = new double *[n];
    for (int i = 0; i < n; ++i) {
        dist[i] = new double[n];
    }
    for (int i = 0; i < nlines; ++i) {
        dist[link[i].u][link[i].v] = dist[link[i].v][link[i].u] = link[i].d;
    }
}

void calcRho()
{
    double percent = 2;
    printf("average percentage of neighbours (hard coded): %5.6lf\n", percent);
    
    int possition = round(nlines * percent / 100);
    sort(link, link + nlines, cmpLink);
    dc = link[possition].d;
    printf("Computing Rho with gaussian kernel of radius: %12.6lf\n", dc);
    
    for (int i = 0; i < n - 1; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double inc = exp(-pow(dist[i][j] / dc, 2));
            rho[i] += inc;
            rho[j] += inc;
        }
    }
//    for (int i = 0; i < nlines; ++i) {
//        if (link[i].d >= dc) {
//            break;
//        }
//        ++rho[link[i].u];
//        ++rho[link[i].v];
//    }
}

void calcDelta()
{
    double maxd = link[nlines - 1].d;
    
    rhoIdx = new int[n];
    for (int i = 0; i < n; ++i) {
        rhoIdx[i] = i;
    }
    sort(rhoIdx, rhoIdx + n, cmpRhoIdx);
    
    const int &idxMaxRho = rhoIdx[n - 1];
    delta[idxMaxRho] = 0;
    neigh[idxMaxRho] = -1;
    for (int i = 0; i < n - 1; ++i) {
        const int &idxi = rhoIdx[i];
        delta[idxi] = maxd;
        neigh[idxi] = -1;
        for (int j = i + 1; j < n; ++j) {
            if (dist[idxi][rhoIdx[j]] <= delta[idxi]) {
                delta[idxi] = dist[idxi][rhoIdx[j]];
                neigh[idxi] = rhoIdx[j];
            }
        }
        
        delta[idxMaxRho] = max(delta[idxi], delta[idxMaxRho]);
    }
    
    printf("Generated file:DECISION GRAPH\n");
    printf("column 1:Density\n");
    printf("column 2:Delta\n");
    FILE *fp = fopen("DECISION_GRAPH", "w");
    for (int i = 0; i < n; ++i) {
        fprintf(fp, "%6.2f %6.2f\n", rho[i], delta[i]);
    }
    fclose(fp);
}

void cluster()
{
    printf("minimum of rho? ");
    double rhomin;
    scanf("%lf", &rhomin);
    printf("minimum of delta? ");
    double deltamin;
    scanf("%lf", &deltamin);
    
    ncluster = 0;
    for (int i = 0; i < n; ++i) {
        if ((rho[i] > rhomin) && (delta[i] > deltamin)) {
            cl[i] = ncluster++;
        } else {
            cl[i] = -1;
        }
    }
    printf("NUMBER OF CLUSTERS: %d \n", ncluster);
    
    icl = new int[ncluster];
    for (int i = 0; i < n; ++i) {
        if (cl[i] != -1) {
            icl[cl[i]] = i;
        }
    }
    
    printf("Performing assignation\n");
    for (int i = n - 1; i >= 0; --i) {
        if (cl[rhoIdx[i]] == -1) {
            cl[rhoIdx[i]] = cl[neigh[rhoIdx[i]]];
        }
    }
    
    for (int i = 0; i < n; ++i) {
        halo[i] = cl[i];
    }
    double *bordRho = new double[ncluster];
    for (int i = 0; i < ncluster; ++i) {
        bordRho[i] = 0;
    }
    
    for (int i = 0; i < n - 1; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if ((cl[i] != cl[j]) && (dist[i][j] <= dc)) {
                double rhoAver = (rho[i] + rho[j]) / 2;
                bordRho[cl[i]] = max(rhoAver, bordRho[cl[i]]);
                bordRho[cl[j]] = max(rhoAver, bordRho[cl[j]]);
            }
        }
    }
    
    for (int i = 0; i < n; ++i) {
        if (rho[i] < bordRho[cl[i]]) {
            halo[i] = -1;
        }
    }
    delete[] bordRho;
    
    for (int i = 0; i < ncluster; ++i) {
        int nc = 0;
        int nh = 0;
        for (int j = 0; j < n; ++j) {
            if (cl[j] == i) {
                ++nc;
            }
            if (halo[j] == i) {
                ++nh;
            }
        }
        printf("CLUSTER: %d CENTER: %d ELEMENTS: %d CORE: %d HALO: %d \n", i + 1, icl[i] + 1, nc, nh, nc - nh);
    }
    
    FILE *fp = fopen("CLUSTER_ASSIGNATION", "w");
    printf("Generated file:CLUSTER_ASSIGNATION\n");
    printf("column 1:element id\n");
    printf("column 2:cluster assignation without halo control\n");
    printf("column 3:cluster assignation with halo control\n");
    for (int i = 0; i < n; ++i) {
        fprintf(fp, "%d %d %d\n", i + 1, cl[i] + 1, halo[i] + 1);
    }
    fclose(fp);
}

int main()
{
    init();
    
    rho = new double[n];
    delta = new double[n];
    neigh = new int[n];
    
    cl = new int[n];
    halo = new int[n];
    
    calcRho();
    calcDelta();
    cluster();
    
    delete rho;
    delete delta;
    delete neigh;
    
    delete cl;
    delete halo;
    
    delete icl;
    for (int i = 0; i < n; ++i) {
        delete dist[i];
    }
    delete dist;
    
    return 0;
}
