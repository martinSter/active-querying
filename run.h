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
#include <assert.h>

// minimum size of outbreak for ground truth
#define MINOUTBREAK 10

// factor to determine cutoff value for T: [0, NUMT * T]
// hence, from T we go back (NUMT - 1) * T
#define NUMT 2

// Number of Monte-Carlo simulations to be performed from each source
#define NSIM 10000

// Number of experiments
#define NEXP 1

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
	double mu;
    double beta;
    // time until first observation
    unsigned int T;
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
    // INFERENCE TIMES
    unsigned int ntinf, *tinf;
	// HELPERS
    unsigned int nps, *ps; // number of possible sources, array with possible sources (all nodes with at least one out-edge)
	unsigned int ns, *s;
    unsigned int true_source, true_source_idx, true_start, ns_gt, *ground_truth;
    unsigned int fin; // first infected node (fin)
    unsigned int nsources; // number of possible sources
    double rank; // rank of true source
    double inv; // log-likelihood that corresponds to the invalid sources (used in posterior ())
    double maxlog; // max-log-likelihood (used in posterior ())
    double os, r0; // average outbreak size and R0
} GLOBALS;

// struct to define nodes
typedef struct NODE {
	unsigned int deg, *nb; // degree, neighbors
	unsigned int *nc, **t; // ordered number of / list of contact times for bisection search
    unsigned int heap, time; // time is 1st the time of infection (for sorting the heap), then the time of recovery (to check if the node is I or R)
    unsigned int e; // used to indicate whether node is in evid or not (onehop uses this information)
    unsigned int inf_by, inf_time, rec_time; // store infection and recovery time of every node
    // unsigned int *ground_truth; // array over inference times (0 = susceptible, 1 = infectious, 2 = recovered)
    unsigned int true_state; // 0 = susceptible, 1 = infectious, 2 = recovered
    unsigned int map_state; // 0 = susceptible, 1 = infectious, 2 = recovered
    unsigned int visited; // used for deterministic spread
    SIR weighted_pred; // store the final prediction (weighted with posterior)
} NODE;

// struct to define node observations for evidence
typedef struct EVID {
    unsigned int n_evid; // number of nodes in evidence
    unsigned int *nodes; // array of nodes in evidence
} EVID;

// struct to define possible sources
typedef struct INF {
    unsigned int n_sources; // number of possible sources, number of active sources
    unsigned int *sources; // array of possible source nodes
    double **log_liks; // array of arrays with log-likelihoods: [[t1, t2, ...],[t1, t2, ...],...] -> outer array are sources
    double *marginal_sources; // array for marginal posterior probabilities of sources
    double *marginal_T; // array for marginal posterior probabilities of T's
} INF;

// struct to define prior on T
// no need for an array for possible values of T as they are defined by [0, n_durs - 1]
typedef struct PRIOR {
    unsigned int n_durs; // number of possible values for T (defines implicitely the cutoff point)
    double *priors; // prior probabilities on possible values for T
} PRIOR;

// struct to store results from FULL
typedef struct OUTF {
    unsigned int *nactives;
    unsigned int *nsources;
    unsigned int *true_sources;
    unsigned int *ranks;
    double *posteriors;
    unsigned int *css; // credible set sizes
    double *sim_time; // CPU run time of simulation process
} OUTF;

// struct to store results from QUERYING
typedef struct OUTQ {
    unsigned int **ranks;
    unsigned int **nsources;
    unsigned int **css;
} OUTQ;

// struct to store CPU run times
typedef struct OUTT {
    double *run_time;
} OUTT;

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
extern unsigned int exptime ();
extern void read_data (FILE *, unsigned int dir);
extern void find_possible_sources ();
extern void init_prior (PRIOR *pr, unsigned int type);
extern void set_up (OUTQ *out);
extern void clean_up (OUTQ *out);
extern void draw_first_node ();
extern void initialize_evidence (EVID *evid);
extern void update_evidence (unsigned int v, EVID *evid);
extern unsigned int find_max (double array[]);
extern void eval_full (unsigned int nexp, INF *inf, OUTF *outf);
extern void eval_query (unsigned int nexp, unsigned int nquery, INF *inf, OUTQ *o);
extern void export_posterior (unsigned int nexp, unsigned int nq, INF *inf, PRIOR *pr);
extern void export_output (OUTF *outf);
extern void export_output_querying (OUTQ *o, char s[]);
extern void export_run_times (OUTT *ot_maxp, OUTT *ot_ucty, OUTT *ot_akld, OUTT *ot_onehop, OUTT *ot_randq);

// quick.c
extern void quick (unsigned int);

// pcg_rnd.c
extern uint16_t pcg_16 ();
extern uint32_t pcg_32 ();
extern uint32_t pcg_32_bounded ();

// sim.c
extern void ground_truth ();
extern void simulate (unsigned int nexp, INF *inf, PRIOR *pr, OUTF *outf);

// detsir.c
extern void find_active_sources (INF *inf, unsigned int start);

// inf.c
extern void full_inference (unsigned int nexp, INF *inf, PRIOR *pr, OUTF *outf);
extern void infer_query (unsigned int nexp, INF *inf, PRIOR *pr, unsigned int (*fn)(INF *inf, PRIOR *pr, EVID *evid), OUTQ *o, OUTT *ot, unsigned int print);

// query.c
extern unsigned int maxp (INF *inf, PRIOR *pr, EVID *evid);
extern unsigned int ucty (INF *inf, PRIOR *pr, EVID *evid);
extern unsigned int akld (INF *inf, PRIOR *pr, EVID *evid);
extern unsigned int onehop (INF *inf, PRIOR *pr, EVID *evid);
extern unsigned int randq (INF *inf, PRIOR *pr, EVID *evid);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
