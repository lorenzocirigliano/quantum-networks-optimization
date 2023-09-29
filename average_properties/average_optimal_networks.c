#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <limits.h>
#include <unistd.h>
#include <float.h>

#define D 2 //The spatial dimension of the hypercube is set to 2 (i.e. points are generated in a 2D square)
#define eta_0 1. //Set the parameter eta_0 (which appears in the definition of capacitance) to 1, as we have done in the simulations
#define lambda_0 1. //Set lambda_0 to 1 and measure L in units of lambda_0

//GLOBAL VARIABLES

unsigned int N;
double L;
double p;
double alpha, alpha_MIN, alpha_MAX, Delta_alpha;
unsigned int M;
unsigned int BOOLE;

double s, q, e, c;

double *pos;
double **W;
struct weight_label ***Q;
int **nn;
unsigned int **a;

double *capacitance;
double *security;
double *efficience;
double *connectivity;

double **hist_deg;
int *bet;
double **hist_bet;

int *dist;
int *prev;

struct weight_label {   // Structure declaration
  double weight;          
  int x;
  int y;       
};

//FUNCTIONS

void FLAG(){    //Helps in debugging
    fprintf(stderr, "\nFLAG\n");
}

double distance(int A, int B){  //Computes the Euclidean distance between two points
    int i;
    double X = 0;
    for(i = 0; i < D; i++){
        X += pow(pos[D*A+i]-pos[D*B+i], 2.);
    }
    return sqrt(X);
}

void create_graph(){    //Creates the FC weighted network
    int i, j;
    for(i = 0; i < N; i++){
        for(j = 0; j < D; j++){
            pos[D*i+j] = -L/2+L*drand48();  //creates coordinates in the D-dimensional square [-L/2,L/2]^D
        }
    }
    for(i=0; i<N; i++){
        for(j=0; j<i; j++){ 
            W[i][j] = -log(1.-eta_0*exp(-distance(i,j)/lambda_0))/log(2);
            W[j][i] = W[i][j];
        }
    }
    for(i = 0; i < N; i ++){    //Set self-capacitances to infinity
        W[i][i] = DBL_MAX; 
    }
}

double min(double a, double b){ //Returns the minimum between two numbers
    if(a < b){ 
        return a;
    }else{
        return b;
    }
}

int pollack(){  //Implements the Pollack algorithm
    int i, j, k, m = 0, check = 1;
    struct weight_label q_temp;
    q_temp.weight = 0;
    for(i = 0; i < N; i++){
        for(j = 0; j < N; j++){
            Q[0][i][j].weight = Q[0][j][i].weight = W[i][j];
            Q[0][i][j].x = i;
            Q[0][i][j].y = j;
        }
    }
    while(check>0){
        check = 0;
        m++;
        Q = (struct weight_label ***)realloc(Q, (1+m)*sizeof(struct weight_label **));
        Q[m] = (struct weight_label **)calloc(N, sizeof(struct weight_label *));
        for(i = 0; i < N; i++) Q[m][i] = (struct weight_label *)calloc(N, sizeof(struct weight_label));
        for(i = 0 ; i < N; i ++){
            for(j = 0 ; j < N; j ++){
                q_temp.weight = 0;
                for(k = 0; k < N; k ++){
                    if(min(Q[0][i][k].weight, Q[m-1][k][j].weight) >= q_temp.weight){
                        if(Q[0][i][k].weight<Q[m-1][k][j].weight){
                            q_temp = Q[0][i][k];
                        }else{
                            q_temp = Q[m-1][k][j];
                        }
                    }
                }
                Q[m][i][j] = q_temp;
                if (Q[m][i][j].weight-Q[m-1][i][j].weight!=0) check++;
            }
        }
    }
    return m;
}


int get_min(int *unv, int U){   //Finds the node with minimum distance in the unvisited array (see dijkstra)
    int i, temp;
    int temp_dist = INT_MAX;
    for(i = 0; i < U; i++){
        if(dist[unv[i]] <= temp_dist){
            temp_dist = dist[unv[i]];
            temp = i;  
        }
    }
    return temp;
}

void remove_node(int **arr, int *size, int index) { //Remove a node from an array and decrease array size
    int i;

    // shift nodes to the left
    for (i = index; i < *size - 1; i++) {
        (*arr)[i] = (*arr)[i + 1];
    }
    (*size)--;

    // resize array
    if(*size == 0){
        free(*arr);
        *arr = NULL;
        *size = 0;
    } else {
    int *temp = realloc(*arr, (*size) * sizeof(int));
        if (temp == NULL) {
            printf("Error: Could not reallocate memory\n");
            return;
        }
    *arr = temp;
    }
}

void dijkstra(int s){   //Implements Dijkstra algorithm to find topological shortest-paths
  
    int i, j, cur_node, n, U = N;
    int index;
    int dist_cur_node, t_dist;
    int *unv;
    int *marked;
    unv = malloc(N * sizeof(int));
    marked = calloc(N, sizeof(int));

    for(i=0; i < N; i++){
        unv[i] = i;
    }

    for (i=0; i<N; i++){
        prev[i] = -1;
        dist[i] = INT_MAX;
    }
    dist[s] = 0;

    while (U != 0){                   
        index = get_min(unv, U);
        cur_node = unv[index];
        remove_node(&unv, &U, index);
        marked[cur_node] = 1;
        dist_cur_node = dist[cur_node];
        if (dist_cur_node == INT_MAX){ 
            free(unv);
            free(marked);
            return;
        }
        for(j = 1; j <= nn[cur_node][0]; j++){
            n = nn[cur_node][j];
            if(marked[n] == 0){
                t_dist = dist_cur_node + 1;
                if (t_dist < dist[n]){
                    dist[n] = t_dist; 
                    prev[n] = cur_node;
                }
            }
        }
    }
    free(marked);
    free(unv);
    return;
}

void build_optimal_path(int start, int end, int x, int y, double q_star){   //Create the optimal path given ``start'' and ``end'' nodes
    int i, j, end_node;

    nn = (int **)calloc(N, sizeof(int *));
    for(i=0; i<N; i++){
        nn[i]=(int *)calloc(1, sizeof(int));
    }
    for(i = 0; i < N; i++){
        for(j = 0; j < i; j++){
            if(Q[0][i][j].weight >= q_star){
                nn[i][0]++;
                nn[i]=(int *)realloc(nn[i], (nn[i][0]+1)*sizeof(int));
                nn[i][nn[i][0]]=j;
                nn[j][0]++;
                nn[j]=(int *)realloc(nn[j], (nn[j][0]+1)*sizeof(int));
                nn[j][nn[j][0]]=i;
            }
        }
    }

    a[x][y] = a[y][x] = 1;

    dist = (int *)calloc(N, sizeof(int));
    prev = (int *)calloc(N, sizeof(int));

    dijkstra(start);
    if(dist[x] < dist[y]){
        end_node = x;
    }else{
        end_node = y;
    }
    i = end_node;
    while(prev[i]!=-1){
        a[i][prev[i]] = a[prev[i]][i] =  1;
        i = prev[i];
        if(i != start && i != end) bet[i]++;
    }

    dijkstra(end);
    if(end_node == x){
        end_node = y;
    }else{
        end_node = x;
    }
    i = end_node;

    while(prev[i]!=-1){
        a[i][prev[i]] = a[prev[i]][i] =  1;
        i = prev[i];
        bet[i]++;
        if(i != start && i != end) bet[i]++;
    }

    for(j=0; j<N; j++){
        free(nn[j]);
    }
    free(nn);

    free(prev);
    free(dist);
}

void pollack_star(int m_max, double alpha, int tau){    //Implements the modified Pollack algorithm
    int i, j, k, m_star, x, y, degree;
    double q_star, e_star, e_check = 0;
    bet = (int *)calloc(N, sizeof(int));
    for(i = 0; i < N; i++){
        for(j = 0; j < i; j++){
            e_star = -DBL_MAX;
            for(k = 0; k < m_max; k++){
                e_check = (1-alpha)*Q[k][i][j].weight + alpha*log(1-p)*k;
                if(e_check > e_star){
                    q_star = Q[k][i][j].weight;
                    x = Q[k][i][j].x;
                    y = Q[k][i][j].y;
                    e_star = e_check;
                    m_star = k;
                }
            }
            if(BOOLE == 1){
                build_optimal_path(i, j, x, y, q_star);
            }
            q += 2*q_star/(N*(N-1));
            s += -2*log(1-p)*m_star/(N*(N-1));
            e += (1-alpha)*2*q_star/(N*(N-1)) + 2*alpha*log(1-p)*m_star/(N*(N-1));
        }
    }
    if(BOOLE==1){
        for(i = 0; i < N; i++){
            degree = 0;
            for(j = 0; j < N; j++){
                c+=(1.*a[i][j])/(N*(N-1));
                degree += a[i][j];
            }
            hist_deg[tau][degree] += 1./(N*M);
            hist_bet[tau][bet[i]] += 1./M;
        }
        free(bet);
    }
}

//MAIN

int main(int argc, char *argv[]){

    if(argc!=9){
        fprintf(stderr, "\nInvalid number of arguments in %s.\nSintax must be: number_of_trusted_nodes\t\tL/lambda_0\t\tp\t\talpha_MIN\t\talpha_MAX\t\tDelta_alpha\t\tNumber_of_realizations\t\tBOOLE\n\n", argv[0]);
        exit(1);
    }

    //Takes input values
    char *useless;
    N = strtol(argv[1], &useless, 10);
    L = strtod(argv[2], &useless);
    p = strtod(argv[3], &useless);
    alpha_MIN = strtod(argv[4], &useless);
    alpha_MAX = strtod(argv[5], &useless);
    Delta_alpha = strtod(argv[6], &useless);    
    M = strtol(argv[7], &useless, 10);
    BOOLE = strtol(argv[8], &useless, 10);
    fprintf(stderr, "\nRunning program with: \nN=%d, \nL/lambda_0=%lf, \np=%lf, \nalpha_MIN=%lf, \nalpha_MAX=%lf, \nDelta_alpha=%lf, \nM=%d.\n\n", N, L, p, alpha_MIN, alpha_MAX, Delta_alpha, M);
    if(N<=0||L<0||p<0||p>1||alpha_MIN<0||alpha_MIN>1||alpha_MAX<0||alpha_MAX>1||Delta_alpha<0||Delta_alpha>1||M<=0||BOOLE<0||BOOLE>1){
        fprintf(stderr, "ERROR: invalid parameters. Check that: \n(1) N>0 (the number of trusted nodes must be a positive integer)\n(2) L>0 (L is a length, measured in units of lambda_0),\n(3) 0<=p<=1 (p is a probability),\n(4) 0<=alpha_MIN,alpha_MAX,Delta_alpha<=1 (alpha ranges at most between 0 and 1),\n(5) M>=1, (the number of realizations must be a positive integer\n(6) BOOLE must be either 0 or 1.\n\n");
        exit(1);
    }

    if(BOOLE==0){
        fprintf(stderr, "REMARK: you decided to NOT COMPUTE topological properties of the optimal network, such as link density, degree distribution and betweenness. This option speeds up the program.\nPlease re-execute the program setting BOOLE=1 otherwise.\n\n");
    }else{
        fprintf(stderr, "REMARK: you decided to COMPUTE the topological properties of the optimal network, such as link density, degree distribution and betweenness. This option slows down the program.\nPlease re-execute the program setting BOOLE=0 otherwise.\n\n");
    }


    //Variables declaration
    int i, j, t, T, m_max, tau;
    double alpha;
    FILE *f1, *f2, *f3;
    
    char file_name_1[50];
    char file_name_2[50];
    char file_name_3[50];


    //Memory allocations
    pos = (double*)calloc(D*N, sizeof(double));
    T = (int)1./Delta_alpha+1;
    capacitance = (double*)calloc(T, sizeof(double));
    security = (double*)calloc(T, sizeof(double));
    efficience = (double*)calloc(T, sizeof(double));
    connectivity = (double*)calloc(T, sizeof(double));
    
    if(BOOLE == 1) {
        hist_deg =  (double **)calloc(T, sizeof(double *));
        hist_bet = (double **)calloc(T, sizeof(double *));
        for(i = 0; i < T; i++){
            hist_deg[i] = (double *)calloc(N, sizeof(double));
            hist_bet[i] = (double *)calloc(N*N, sizeof(double));
        }
    }
    
    W = (double **)calloc(N, sizeof(double *));
    for(i = 0; i < N; i++){
        W[i] = (double *)calloc(N, sizeof(double));
    }

    srand48(time(NULL));    //Initialize RNG with time(NULL)
	for(i=0; i<1000; i++){
        drand48();
	}

    //Define output files names
    sprintf(file_name_1, "efficiency_N=%d_L=%.4lf_p=%.2lf.dat", N, L/lambda_0, p);
    if((f1 = fopen(file_name_1, "w")) == NULL) {
        fprintf(stderr, "Could not open a file. Process interrupted.\n");
        exit(1);
    } 

    if(BOOLE==1){
        sprintf(file_name_2, "degree_distribution_N=%d_L=%.4lf_p=%.2lf.dat", N, L/lambda_0, p);
        sprintf(file_name_3, "betweenness_distribution_N=%d_L=%.4lf_p=%.2lf.dat", N, L/lambda_0,p);
        if((f2 = fopen(file_name_2, "w")) == NULL) {
            fprintf(stderr, "Could not open a file. Process interrupted.\n");
            exit(1);
        }  
        if((f3 = fopen(file_name_3, "w")) == NULL) {
            fprintf(stderr, "Could not open a file. Process interrupted.\n");
            exit(1);
        } 
    }

    fprintf(f1, "#alpha\t\tQ\t\t\tS\t\t\tE\t\t\tC\n");
    e = s = q = c = 0;
    for(t = 1; t <= M; t++){ //Run the optimization algorithm for M independent realizations of points in the square
        Q = (struct weight_label ***)calloc(1, sizeof(struct weight_label **));
        Q[0] = (struct weight_label **)calloc(N, sizeof(struct weight_label *));
        for(i=0; i<N; i++){
            Q[0][i] = (struct weight_label *)calloc(N, sizeof(struct weight_label));
        }
        create_graph();
        m_max = pollack();
        tau = 0;
        for(alpha = alpha_MIN; alpha <= alpha_MAX; alpha += Delta_alpha){   //Find average optimal network's properties for each alpha
            if(BOOLE == 1){
                a = (unsigned int **)calloc(N, sizeof(unsigned int *));
                for(i=0; i<N; i++){
                    a[i]=(unsigned int *)calloc(N, sizeof(unsigned int));
                }
            }
            pollack_star(m_max, alpha, tau);
            capacitance[tau] += q/M;
            security[tau] += s/M;
            efficience[tau] += e/M;
            connectivity[tau] += c/M;
            if(BOOLE == 1){
                for(i=0; i<N; i++){
                    free(a[i]);
                }
                free(a);
            }
            e = s = q = c = 0;
            tau++;
        }
        for(i=0; i<=m_max; i++){
            for(j=0; j<N; j++){
                free(Q[i][j]);
            }
            free(Q[i]);
        }
        free(Q);
    }
    tau = 0;
    
    for(alpha = alpha_MIN; alpha <= alpha_MAX; alpha += Delta_alpha){
        fprintf(f1, "%lf\t%lf\t%lf\t%lf\t%lf\n", alpha, capacitance[tau], security[tau], efficience[tau], connectivity[tau]);
        if(BOOLE == 1){
            fprintf(f2, "#alpha=%.2lf\n", alpha);
            fprintf(f3, "#alpha=%.2lf\n", alpha);
            for(i = 0; i < N; i++){
                fprintf(f2, "%d\t%lf\n", i, hist_deg[tau][i]);
            }
            for(i = 0; i < N*N; i++){
                fprintf(f3, "%d\t%lf\n", i, hist_bet[tau][i]);
            }
        }
        tau++;
    }

    //Free allocated memory
    for(i=0; i<N; i++){
        free(W[i]);
    }
    free(W);

    free(pos);
    free(capacitance);
    free(security);
    free(efficience);
    free(connectivity);
        
    fclose(f1);

    if(BOOLE == 1){
        for(i=0; i<T; i++){
            free(hist_deg[i]);
            free(hist_bet[i]);
        }
        free(hist_deg);
        free(hist_bet);
        fclose(f2);
        fclose(f3);        
    }
    return(0);
}
