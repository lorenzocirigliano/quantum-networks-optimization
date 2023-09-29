#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <limits.h>
#include <unistd.h>
#include <float.h>
#include <sys/types.h>
#include <sys/wait.h>

#define D 2 //The spatial dimension of the hypercube is set to 2 (i.e. points are generated in a 2D square)
#define eta_0 1.    //Set the parameter eta_0 (which appears in the definition of capacitance) to 1, as we have done in the simulations
#define alpha_MIN 0.    //Set the minimum value of alpha for optimization to 0
#define alpha_MAX 1.    //Set the largest value of alpha for optimization to 1

//GLOBAL VARIABLES

unsigned int N;
double lambda_0;
double p;
double alpha, Delta_alpha;

char input_file_name[50];
char output_file_name[150];

double s, q, e, c;

double *pos;
double **W;
struct weight_label ***Q;
double **E;
int **nn;
unsigned int **a;

double *capacitance;
double *security;
double *efficience;

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

void create_graph(){
    int i, j;
    for(i=0; i<N; i++){
        for(j=0; j<i; j++){ //Creates the FC weighted network
            W[i][j] = -log(1.-exp(-distance(i,j)/lambda_0))/log(2);
            W[j][i] = W[i][j];
        }
    }
    for(i = 0; i < N; i ++){    //Set self-capacitances to infinity
        W[i][i] = DBL_MAX; 
    }
}

//This function reads as input the coordinates (x,y) of N points in a square  
void read_input(){
    int i;
    double x, y;
    FILE *f1;

    if((f1 = fopen(input_file_name, "r")) == NULL) {
        fprintf(stderr, "Could not open a file. Process interrupted.\n");
        exit(1);
    }  
    i = 0;
    while(fscanf(f1, "%lf %lf", &x, &y) == 2){
        pos[2*i] = x;
        pos[2*i+1] = y;
        i++;
    }
    fclose(f1);
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
    }

    for(j=0; j<N; j++){
        free(nn[j]);
    }
    free(nn);

    free(prev);
    free(dist);
}

void pollack_star(int m_max, int BOOLE){    //Implements the modified Pollack algorithm
    int i, j, k, m_star, x, y;
    double q_star, e_star, e_check = 0;
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
            if(BOOLE == 1) build_optimal_path(i, j, x, y, q_star);
            q += 2*q_star/(N*(N-1));
            s += -2*log(1-p)*m_star/(N*(N-1));
            e += (1-alpha)*2*q_star/(N*(N-1)) + 2*alpha*log(1-p)*m_star/(N*(N-1));
        }
    }
}

void print_graph(){
    int i, j;
    FILE *f1;

    sprintf(output_file_name, "list_of_connections_lambda_0=%.4lf_p=%.2lf_alpha=%lf_for_%s", lambda_0, p, alpha, input_file_name);
    if((f1 = fopen(output_file_name, "w")) == NULL) {
        fprintf(stderr, "Could not open a file. Process interrupted.\n");
        exit(1);
    } 

    for(i = 0; i < N; i++){
        for(j = 0; j < i; j++){
            if(a[i][j] == 1) fprintf(f1, "%d %d\n", i, j);
        }
    }
    fclose(f1);
}

int main(int argc, char *argv[]){

    if(argc!=6) {
        fprintf(stderr, "\nInvalid number of arguments in %s.\nSintax must be: <input_file_name>\t\tnumber_of_trusted_nodes\t\tlambda_0\t\tp\t\tDelta_alpha\n\n", argv[0]);
        exit(1);
    }
    //Takes input values
    char *useless;
    strcpy(input_file_name, argv[1]);    
    N = strtol(argv[2], &useless, 10);
    lambda_0 = strtod(argv[3], &useless);
    p = strtod(argv[4], &useless);
    Delta_alpha = strtod(argv[5], &useless);  
   
    //PART 1
    fprintf(stderr, "\n(1) Running part 1 with N=%d, lambda_0=%.2lf, p=%lf,Delta_alpha=%lf\n\n", N, lambda_0, p, Delta_alpha);

    if(N<=0||lambda_0<0||p<0||p>1||Delta_alpha<0||Delta_alpha>1){
        fprintf(stderr, "ERROR: invalid parameters. Check that: \n(1) N>0 (the number of trusted nodes must be a positive integer)\n(2) lambda_0>0 (lambda_0 is a length),\n(3) 0<=p<=1 (p is a probability),\n(4) 0<=Delta_alpha<=1 (alpha ranges at most between 0 and 1).\n\n");
        exit(1);
    }    

    //Variables declaration
    int i, j, T, m_max, tau;
    FILE *f1;

    //Memory allocations
    pos = (double*)calloc(D*N, sizeof(double));
    T = (int)1./Delta_alpha+1;
    capacitance = (double*)calloc(T, sizeof(double));
    security = (double*)calloc(T, sizeof(double));
    efficience = (double*)calloc(T, sizeof(double));

    W = (double **)calloc(N, sizeof(double *));
    E = (double **)calloc(N, sizeof(double **));
    for(i = 0; i < N; i++){
        W[i] = (double *)calloc(N, sizeof(double));
        E[i] = (double *)calloc(N, sizeof(double));
    }

    read_input(); 
    create_graph();

    sprintf(output_file_name, "optimal_efficiency_for_%s", input_file_name);
    if((f1 = fopen(output_file_name, "w")) == NULL) {
        fprintf(stderr, "Could not open a file. Process interrupted.\n");
        exit(1);
    }  
    //fprintf(f1, "#alpha\t\tQ\t\t\tL\t\t\tS\t\t\tE\n");
    fprintf(f1, "#alpha\t\tQ*\t\t\tL*\n");
    e = s = q = c = 0;

    Q = (struct weight_label ***)calloc(1, sizeof(struct weight_label **));
    Q[0] = (struct weight_label **)calloc(N, sizeof(struct weight_label *));
    for(i=0; i<N; i++){
        Q[0][i] = (struct weight_label *)calloc(N, sizeof(struct weight_label));
    }
    m_max = pollack();
    tau = 0;
    for(alpha = alpha_MIN; alpha <= alpha_MAX; alpha += Delta_alpha){   //Find average optimal network's properties
        pollack_star(m_max, 0); //Here the second entrance is set to 0 since we want only to find efficience and not the optimal paths
        capacitance[tau] += q;
        security[tau] += s;
        efficience[tau] += e;
        e = s = q = 0;
        tau++;
    }
    for(i=0; i<=m_max; i++){
        for(j=0; j<N; j++){
            free(Q[i][j]);
        }
        free(Q[i]);
    }
    free(Q);
    tau = 0;

    for(alpha = alpha_MIN; alpha <= alpha_MAX; alpha += Delta_alpha){
        //fprintf(f1, "%lf\t%lf\t%lf\t%lf\t%lf\n", alpha, capacitance[tau], -security[tau]/log(1-p)+1, security[tau], efficience[tau]);
        fprintf(f1, "%lf\t%lf\t%lf\n", alpha, capacitance[tau], -security[tau]/log(1-p)+1);
        tau++;
    }
    
    free(capacitance);
    free(security);
    free(efficience);
    fclose(f1);
    fprintf(stderr, "Part 1 is terminated.\nPlease inspect the file '%s'to decide the alpha value for which you want the optimal network.\n\n***\n\n***\n", output_file_name);

    //PART 2
    fprintf(stderr, "\n(2) Running part 2 of the program.\n\nInsert the chosen value of alpha: ");
    scanf("%lf", &alpha);
    fprintf(stderr, "\n");
    tau = 0;
    while(alpha<0||alpha>1){
        fprintf(stderr, "ERROR: invalid value of alpha. It must be between 0 and 1. Insert again: ");
        scanf("%lf", &alpha);
        fprintf(stderr, "\n");
        tau ++;
        if(tau == 5){
            fprintf(stderr, "Too many errors. Program aborted.\n\n");
            exit(1);
        }
    } 

    e = s = q = c = 0;
    Q = (struct weight_label ***)calloc(1, sizeof(struct weight_label **));
    Q[0] = (struct weight_label **)calloc(N, sizeof(struct weight_label *));
    for(i=0; i<N; i++){
        Q[0][i] = (struct weight_label *)calloc(N, sizeof(struct weight_label));
    }
    create_graph();
    m_max = pollack();

    //Build the network for the desired value of alpha
    a = (unsigned int **)calloc(N, sizeof(unsigned int *));
    for(i=0; i<N; i++){
        a[i]=(unsigned int *)calloc(N, sizeof(unsigned int));
    }

    pollack_star(m_max, 1); //Here the second entrance is set to 1 in order to build the network
    print_graph();

    for(i=0; i<N; i++){
        free(a[i]);
    }
    free(a);  
    for(i = 0; i <= m_max; i++){
        for(j=0; j<N; j++){
            free(Q[i][j]);
        }
        free(Q[i]);
    }
    free(Q);
    fprintf(stderr, "Part 2 is terminated. You can find the list of connections of the optimal network in the file '%s'.\n\n", output_file_name);
    //Free allocated memory
    for(i=0; i<N; i++){
        free(W[i]);
    }
    free(W);

    free(pos);

    return(0);
}
