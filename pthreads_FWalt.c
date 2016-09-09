/*************************************************************************/
// gcc -pthread -o testPthreadsAlt pthreads_FW.c -lm -lrt

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <pthread.h>

#define INFINITE 100000000
#define DISTANCE 10           // determines max weight of the edges  
#define EDGES 1048             // determines density of graph unless graph is complete 
#define VERTICES 8000


#define GIG 1000000000
#define CPG 1.4           // Cycles per GHz -- Adjust to your computer

#define BASE VERTICES
#define ITERS 1
#define DELTA 1

#define OPTIONS 2
#define IDENT 0

#define MINVAL   0.0
#define MAXVAL  DISTANCE

#define TOL 1e-6

typedef int data_t; 


/* Create abstract data type for matrix */
typedef struct {
  long int len;
  long int edges;
  data_t *data;
} matrix_rec, *matrix_ptr;

int NUM_THREADS = 8;

/* used to pass parameters to worker threads */
struct thread_data{
  int thread_id;
  matrix_ptr G;
  matrix_ptr pi;
};

pthread_barrier_t barrier;


/************************************************************************/

int main(int argc, char* argv[]){
  int OPTION;
  struct timespec diff(struct timespec start, struct timespec end);
  struct timespec time1, time2;
  struct timespec time_stamp[OPTIONS][ITERS+1];
  int clock_gettime(clockid_t clk_id, struct timespec *tp);
  matrix_ptr new_matrix(long int len);
  int set_matrix_length(matrix_ptr m, long int index);
  long int get_matrix_length(matrix_ptr m);
  int set_matrix_edges(matrix_ptr G, long int edges);
  long int get_matrix_edges(matrix_ptr G);
  int init_matrix(matrix_ptr m, long int len);
  int zero_matrix(matrix_ptr m, long int len);

  void *cb_work(void *threadarg);
  void pt_fw(matrix_ptr G, matrix_ptr pi);

  void compareResults(matrix_ptr pt_fw, matrix_ptr serial_fw, int vertices);

  //void init_distance_matrix(matrix_ptr G, long int len, long int distance, long int edges, long int maxedges); // use when variable density is required
  // void init_distance_matrix(matrix_ptr G, matrix_ptr G_pt, long int len, long int distance, int seed);
  void init_distance_matrix(matrix_ptr G, long int len, long int distance); // use for complete graphs
  void init_predecessor(matrix_ptr pi, long int len);
  void print_graph(matrix_ptr G);
  void serial_floyd_warshall(matrix_ptr G, matrix_ptr pi);

  long int i, j, k;
  long int time_sec, time_ns;
  long int MAXSIZE = BASE+(ITERS+2)*DELTA;
  long int MAXEDGES = EDGES+(ITERS+2)*DELTA; 

  printf("\nAllocating matrices\n");
  matrix_ptr G = new_matrix(MAXSIZE);
  matrix_ptr pi = new_matrix(MAXSIZE);
  matrix_ptr G_pt = new_matrix(MAXSIZE);
  matrix_ptr pi_pt = new_matrix(MAXSIZE);
  //init_distance_matrix(G, MAXSIZE, (int)DISTANCE, (int)EDGES, MAXEDGES);
  //init_predecessor(pi, MAXSIZE);

  printf("\nEntering options cycle\n");

  // NUM_THREADS = 8;
  // for each iteration, initialize new distance and predecessor matrices and then time the algorithm's execution
  OPTION = 0;
  for (i = 0; i < ITERS; i++) {
    //init_distance_matrix(G,BASE+(i+2)*DELTA, DISTANCE, EDGES+(i+2)*DELTA, MAXEDGES);
    init_distance_matrix(G, BASE+(i+2)*DELTA, DISTANCE); 
    init_predecessor(pi,BASE+(i+2)*DELTA);
    clock_gettime(CLOCK_REALTIME, &time1);
    serial_floyd_warshall(G, pi);
    pt_fw(G,pi);
    clock_gettime(CLOCK_REALTIME, &time2);
    time_stamp[OPTION][i] = diff(time1,time2);
    // printf("\niter = %d", i);
  }

  // printf("\nnext Option\n");

  NUM_THREADS = 8;
  OPTION++;
  for (i = 0; i < ITERS; i++) {
    init_distance_matrix(G,BASE+(i+2)*DELTA, DISTANCE);
    init_predecessor(pi,BASE+(i+2)*DELTA);
    clock_gettime(CLOCK_REALTIME, &time1);
    printf("\nhere\n");
    pt_fw(G_pt,pi);
    clock_gettime(CLOCK_REALTIME, &time2);
    time_stamp[OPTION][i] = diff(time1,time2);
    printf("\niter = %d", i);
  }


  printf("\nvertices, serial_fw, pt_fw");
  for (i = 0; i < ITERS; i++) {
    printf("\n%d, ", BASE+(i+2)*DELTA);
    for (j = 0; j < OPTIONS; j++) {
      if (j != 0) printf(", ");
      printf("%ld", (long int)((double)(CPG)*(double)
		 (GIG * time_stamp[j][i].tv_sec + time_stamp[j][i].tv_nsec)));
    }
  }

  printf("\n\nTime: \n");
  printf("\nvertices, pt_fw");
  for (i = 0; i < ITERS; i++) {
    printf("\n%d, ", BASE+(i+2)*DELTA);
    for (j = 0; j < OPTIONS; j++) {
      if (j != 0) printf(", ");
      printf("%ld", (long int)((double)
		 (time_stamp[j][i].tv_sec + time_stamp[j][i].tv_nsec)));
    }
  }

  
  //print_graph(G);
  //print_graph(G_pt);
  //print_graph(pi);

  printf("\n");

  // compareResults(G_pt, G, VERTICES);

  printf("\n");

  return 0;
}






/**********************************************/

/* Create matrix of specified length */
matrix_ptr new_matrix(long int len)
{
  long int i;

  /* Allocate and declare header structure */
  matrix_ptr result = (matrix_ptr) malloc(sizeof(matrix_rec));
  if (!result) return NULL;  /* Couldn't allocate storage */
  result->len = len;

  /* Allocate and declare array */
  if (len > 0) {
    data_t *data = (data_t *) calloc(len*len, sizeof(data_t));
    if (!data) {
	  free((void *) result);
	  printf("\n COULDN'T ALLOCATE STORAGE \n", result->len);
	  return NULL;  /* Couldn't allocate storage */
	}
	result->data = data;
  }
  else result->data = NULL;

  return result;
}

/* Set length of matrix */
int set_matrix_length(matrix_ptr m, long int index)
{
  m->len = index;
  return 1;
}


/*Set number of edges of graph*/
int set_matrix_edges(matrix_ptr G, long int edges){
  G->len = edges;
  return 1;
}


/* Return length of matrix */
long int get_matrix_length(matrix_ptr m)
{
  return m->len;
}

/* Return number of edges of graph */
long int get_matrix_edges(matrix_ptr G){
  return G->edges;
}

/* initialize matrix */
int init_matrix(matrix_ptr m, long int len)
{
  long int i;

  if (len > 0) {
    m->len = len;
    for (i = 0; i < len*len; i++)
      m->data[i] = (data_t)(i);
    return 1;
  }
  else return 0;
}

/* initialize matrix */
int zero_matrix(matrix_ptr m, long int len)
{
  long int i,j;

  if (len > 0) {
    m->len = len;
    for (i = 0; i < len*len; i++)
      m->data[i] = (data_t)(IDENT);
    return 1;
  }
  else return 0;
}

data_t *get_matrix_start(matrix_ptr m)
{
  return m->data;
}

/*************************************************/

struct timespec diff(struct timespec start, struct timespec end)
{
  struct timespec temp;
  if ((end.tv_nsec-start.tv_nsec)<0) {
    temp.tv_sec = end.tv_sec-start.tv_sec-1;
    temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
  } else {
    temp.tv_sec = end.tv_sec-start.tv_sec;
    temp.tv_nsec = end.tv_nsec-start.tv_nsec;
  }
  return temp;
}

/*************************************************/

double fRand(double fMin, double fMax)
{
    double f = (double)random() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

/*************************************************/

void print_graph(matrix_ptr G)
{
  int i, j;
  int len = get_matrix_length(G);
  data_t* graph = get_matrix_start(G);
  for(i = 0; i < len; i++){
    for(j = 0; j < len; j++) {	
      if (graph[i * len + j] < INFINITE){
	printf("%d ", graph[i * len + j]);
      }
      else{
	printf("INF ");
      }
    }
    printf("\n");
  }
  printf("\n");
}


// generates a 2D distance matrix for a complete directed graph
void init_distance_matrix(matrix_ptr G, long int len, long int distance){
  printf("\nInitialize matrix: %d vertices\n", len);
  srand(time(0));
  long int i, j;
  G->len = len;
  G->edges = len * len;
  for(i = 0; i < len; i++){
    for(j = 0; j < len; j++){
      if(i == j){
	G->data[ i * len + j] = (data_t)0;
	continue;
      }
      long int edgeWeight = rand() % distance;  
      G->data[i*len+j] = edgeWeight == 0 ? 1 : edgeWeight;
    }
  }
  //print_graph(G);
}


/*
// generates a 2D distance matrix for a complete directed graph
void init_distance_matrix(matrix_ptr G, matrix_ptr G_pt, long int len, long int distance){
  printf("\nInitialize matrix: %d vertices\n", len);
  srand(time(0));
  long int i, j;
  G->len = len;
  G_pt->len = len;
  G->edges = len * len;
  G_pt->edges = len * len;
  for(i = 0; i < len; i++){
    for(j = 0; j < len; j++){
      if(i == j){
	G->data[ i * len + j] = (data_t)0;
	G_pt->data[ i * len + j] = (data_t)0;
	continue;
      }
      long int edgeWeight = rand() % distance;  
      G->data[i*len+j] = edgeWeight == 0 ? 1 : edgeWeight;
      G_pt->data[i*len+j] = edgeWeight == 0 ? 1 : edgeWeight;
    }
  }
  //print_graph(G);
}
*/

/*  // uncomment function to intialize an incomplete graph with parameterized density
void init_distance_matrix(matrix_ptr G, long int len, long int distance, long int edges, long int maxedges, int seed){
  srand(seed);
  long int edgeRange = (maxedges/edges); 
  long int i, j;
  G->len = len;
  G->edges = edges;
  for(i = 0; i < len; i++){
    for(j = 0; j < len; j++){
      if(i == j){
	G->data[ i * len + j] = (data_t)0;
	continue;
      }
      long int edgeWeight = edgeRange; 
      G->data[i*len+j] = edgeWeight == 0 ? (rand()%distance)+1)  : INFINITE; //set edge random edge weight
    }
  }
  //print_graph(G);
} 
*/


// generates a 2D predecessor matrix for a directed graph
void init_predecessor(matrix_ptr pi, long int len){
  pi->len = len;
  int i, j;
  for(i = 0; i < len; i++){
    for(j = 0; j < len; j++){
      pi->data[i*len+j] = (data_t)(-1); 
    }
  }

  //print_graph(pi);

}


/*************************************************/

void serial_floyd_warshall(matrix_ptr G, matrix_ptr pi){
  long int vertices = get_matrix_length(G);
  data_t* distanceMatrix = get_matrix_start(G);
  data_t* predMatrix = get_matrix_start(pi);
  long int i, j, k;

  for(k = 0; k < vertices; k++){
    for(i = 0; i < vertices; i++){
      for(j = 0; j < vertices; j++){
	int currentIndex = i*vertices+j;
	int currentPath = distanceMatrix[currentIndex];
	int nextPath = distanceMatrix[i*vertices+k] + distanceMatrix[k*vertices+j];
	if(nextPath < currentPath){
	  distanceMatrix[currentIndex] = nextPath;
	  predMatrix[currentIndex] = k;
	}
      }
    }

    //  print_graph(G);
    //  print_graph(pi);
  }

}


/* CPU bound multithreaded code*/
/* first, the actual thread function */
void *cb_work(void *threadarg)
{
  long int i, j, k, low, high;
  struct thread_data *my_data;
  my_data = (struct thread_data *) threadarg;
  int taskid = my_data->thread_id;
  matrix_ptr G = my_data->G;
  matrix_ptr pi = my_data->pi;
  long int vertices = get_matrix_length(G);
  data_t* distanceMatrix = get_matrix_start(G);
  data_t* predMatrix = get_matrix_start(pi);
  

  low = (taskid * vertices)/NUM_THREADS;
  high = ((taskid+1)* vertices)/NUM_THREADS;
 
  for(k = 0; k < vertices; k++){

    //  printf("\nwaiting for barrier\n");

    pthread_barrier_wait(&barrier);

    for(i = low; i < high; i++){
      for(j = 0; j < vertices; j++){
	int currentIndex = i*vertices+j;
	int currentPath = distanceMatrix[currentIndex];
	int nextPath = distanceMatrix[i*vertices+k] + distanceMatrix[k*vertices+j];
	if(nextPath < currentPath){
	  distanceMatrix[currentIndex] = nextPath;
	  predMatrix[currentIndex] = k;
	}
      }
    }
  }

  pthread_exit(NULL);
}

/* Now, the pthread calling function */
void pt_fw(matrix_ptr G, matrix_ptr pi)
{
  long int i, j, k;
  pthread_t threads[NUM_THREADS];
  struct thread_data thread_data_array[NUM_THREADS];
  pthread_barrier_init(&barrier, NULL, NUM_THREADS);
  int rc;
  long t;

  printf("\ngenerating threads...\n");

  for (t = 0; t < NUM_THREADS; t++) {
    thread_data_array[t].thread_id = t;
    thread_data_array[t].G = G;
    thread_data_array[t].pi = pi;
    rc = pthread_create(&threads[t], NULL, cb_work,
			(void*) &thread_data_array[t]);
    if (rc) {
      printf("ERROR; return code from pthread_create() is %d\n", rc);
      exit(-1);
    }
  }

  for (t = 0; t < NUM_THREADS; t++) {
    if (pthread_join(threads[t],NULL)){ 
      printf("ERROR; code on return from join is %d\n", rc);
      exit(-1);
    }
  }
}

void compareResults(matrix_ptr pt_fw, matrix_ptr serial_fw, int vertices){
  long int i, errCount = 0;
  data_t* pt_dm = get_matrix_start(pt_fw);
  data_t* serial_dm = get_matrix_start(serial_fw);
  // Compare the results
  for(i = 0; i < vertices * vertices; i++) {
    if (abs(pt_dm[i] - serial_dm[i]) > TOL) {
      errCount++;
    }
  }
	
  /*
    for(i = 0; i < 100; i++) {
    printf("%d:\t%.8f\t%.8f\n", i, h_result_gold[i], h_result[i]);
    }
  */
	
  if (errCount > 0) {
    printf("\n@ERROR: TEST FAILED: %d results did not matched\n", errCount);
  }
        


  else {
    printf("\nTEST PASSED: All results matched\n");
  }
}
