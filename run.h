// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// header file

// inspired by code for temporal network SIR by Petter Holme (2018)

// include external libraries
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdint.h>
#include <string.h>
#include <float.h>

// minimum size of outbreak for ground truth
#define MINOUTBREAK 10

// end of time period
#define TEND 5

// time steps for inference times
// first inference is t0 + TSTEPS
// second inference time is t0 + 2*TSTEPS
#define TSTEPS 1

// number of Monte-Carlo simulations to be performed from each source
#define NSIM 10000

// number of experiments
#define NEXP 10

// number of queries
#define NQUERIES 10

// NONE and END are used in various ways, the only purpose of NONE < END is for the S(x) macro
// UINT_MAX is the maximum value for an object of type unsigned int.
// Those integers are effectively used to assign time infinity.
#define NONE (UINT_MAX - 1)
#define END UINT_MAX

// INVALID will be set in log-likelihoods for invalid sources
#define INVALID 100

// is x susceptible?
#define S(x) (n[(x)].heap < END)

// auxiliary macro to compute square of x
#define SQ(x) ((x) * (x))

// define precision
#define PREC pow(10, -16)

// macro to avoid unused parameter warning
#define UNUSED(x) (void)(x)

// struct to define SIR probabilities for one node
typedef struct SIR {
    double s;
    double i;
    double r;
} SIR;

// struct to define global parameters
typedef struct GLOBALS {
	// SIR INPUT PARAMETERS
    double beta;
    // ARRAYS FOR RANDOM NUMBER GENERATION
	unsigned short rnd2inx[0x10001]; // mapping 16-bit random number to index (beta)
    double random_uniform[0x10001];
	// NETWORK SPECS (number of nodes, number of edges, duration, starting point of epidemic)
	unsigned int n, e, dur, t0;
    // TIME OF INFERENCE
    unsigned int t_now;
	// OTHER GLOBALS FOR HEAP
	unsigned int nheap, *heap;
	// FOR RND (random number generation)
	uint64_t state;
	uint32_t rmem;
	unsigned int cutoff_source, cutoff_dur; // to get the probabilities right . .
    float rexp[0x10000];
    // INFERENCE TIMES
    unsigned int ntinf;
    float *tinf;
	// HELPERS
	unsigned int ns, *s;
    unsigned int true_source, ns_gt, *ground_truth;
    unsigned int fin; // first infected node (fin)
    unsigned int nsources; // number of possible sources
    unsigned int qmap; // q^MAP node
    double rank; // rank of true source
    double inv; // log-likelihood that corresponds to the invalid sources (used in posterior ())
    double maxlog; // max-log-likelihood (used in posterior ())
    double os, r0; // average outbreak size and R0
} GLOBALS;

// struct to define nodes
typedef struct NODE {
	unsigned int deg, *nb; // degree, neighbors
    unsigned int heap; // ...
    float time, inf_time, rec_time;
    unsigned int e; // used to indicate whether node is in evid or not (ONE-HOP uses this information)
    unsigned int inf_by;
    unsigned int *ground_truth; // array over inference times (0 = susceptible, 1 = infectious, 2 = recovered)
    unsigned int true_state; // 0 = susceptible, 1 = infectious, 2 = recovered
    unsigned int map_state; // 0 = susceptible, 1 = infectious, 2 = recovered
    SIR weighted_pred; // store the final prediction (weighted with posterior)
} NODE;

// struct to define node observations for evidence
typedef struct EVID {
    unsigned int n_evid; // number of nodes in evidence
    unsigned int *nodes; // array of nodes in evidence
} EVID;

// struct to define possible sources
typedef struct POSTERIOR {
    double *log_liks; // array with log-likelihoods
    double *posterior; // array with posterior probabilities
} POSTERIOR;

// struct to store results from FULL
typedef struct OUTPUT {
    unsigned int *truesource;
    unsigned int **qmap;
    unsigned int **qrand;
    unsigned int **ranks;
    unsigned int **nsources;
    unsigned int **nactive;
} OUTPUT;

// struct to store results from QUERYING
typedef struct OUTPUTQ {
    unsigned int ***qmap;
    unsigned int ***ranks;
    unsigned int ***nsources;
} OUTPUTQ;

// struct for Monte-Carlo probabilities
typedef struct MCS {
    double **inf; // array of arrays for probabilities of being infectious
    double **rec; // array of arrays for probabilities of being recovered
} MCS;

// FUNCTION PROTOTYPES FROM ALL SOURCE FILES
// the 'extern' statement is not absolutely necessary
// (the compiler implicitly assumes it for functions)

// heap.c
extern void up_heap (unsigned int);
extern void del_root ();

// misc.c
extern void read_data (FILE *);
extern void set_up (OUTPUTQ *out);
extern void clean_up (OUTPUTQ *out);
extern void draw_first_node (unsigned int tinf_idx);
extern void initialize_evidence (EVID *evid);
extern void update_evidence (unsigned int v, EVID *evid);
extern unsigned int find_max (double array[]);
extern void sort_rank (POSTERIOR *post);
extern unsigned int random_draw (unsigned int inf_time_idx);
extern void export_output (OUTPUT *out);
extern void export_output_querying (OUTPUTQ *out, char s[]);

// quick.c
extern void quick (unsigned int);

// pcg_rnd.c
extern uint16_t pcg_16 ();
extern uint32_t pcg_32 ();
extern uint32_t pcg_32_bounded ();

// sim.c
extern void ground_truth ();
extern void monte_carlo ();

// inf.c
extern void full_inference (unsigned int nexp, POSTERIOR *post, OUTPUT *o);
extern void run_querying (unsigned int nexp, POSTERIOR *post, OUTPUTQ *o1, OUTPUTQ *o2, OUTPUTQ *o3, OUTPUTQ *o4, OUTPUTQ *o5);

// query.c
extern unsigned int maxp (unsigned int tinf_idx, POSTERIOR *post, EVID *evid);
extern unsigned int ucty (unsigned int tinf_idx, POSTERIOR *post, EVID *evid);
extern unsigned int akld (unsigned int tinf_idx, POSTERIOR *post, EVID *evid);
extern unsigned int onehop (unsigned int tinf_idx, POSTERIOR *post, EVID *evid);
extern unsigned int randq (unsigned int tinf_idx, POSTERIOR *post, EVID *evid);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
