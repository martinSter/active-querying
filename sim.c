// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// simulation functions

// include header file
#include "run.h"

// declare external variables
extern GLOBALS g;
extern NODE *n;
extern MCS *p;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// this routine first localizes the first contact later than 'now' in t
// then picks the contact that can infect (a chain of bernouilli trials)
// among the rest of the contacts. It returns the time of the infecting contact

unsigned int contagious_contact (unsigned int *t, unsigned int nt, unsigned int now) {
	
    // declare unsigned integers
    unsigned int lo = 0, mid, hi = nt - 1;

    // no need to search further bcoz t is sorted. Note that the bisection search depends on this line.
	if (t[hi] <= now) return END;

	// the actual bisection search
	while (lo < hi) {
		mid = (lo + hi) >> 1;
		if (t[mid] > now) hi = mid;
		else lo = mid + 1;
	}

	// get a random contact
	hi += g.rnd2inx[pcg_16()];

    // if the contact is too late, skip it
	if (hi >= nt) return NONE;

	// return the time of the contact
	return t[hi];
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// this routine does the book keeping for an infection event

void infect () {
    
    // declare unsigned integers
	unsigned int i, you, t, me = g.heap[1];
    
    // set now to infection time of me and sample time until recovery
	unsigned int now = n[me].time, duration = exptime();
    
    // store infection and recovery time of me
    n[me].inf_time = n[me].time;
    n[me].rec_time = n[me].time + duration;

    // take the newly infected off the heap
	del_root();

    // if the duration is zero, no one else can be infected
	if (duration > 0) {
        
        // set time of me to time of recovery
		n[me].time += duration;

		// go through the neighbors of the infected node
		for (i = 0; i < n[me].deg; i++) {
            // set you to neighbor at index i
			you = n[me].nb[i];
            // if you is S, you can be infected
			if (S(you)) {
				// find the infection time of you
				t = contagious_contact(n[me].t[i], n[me].nc[i], now);
				// bcoz of the sorting of nbs, we can break
                if (t == END) break;

				// if the infection time is before when me gets recovered,
				// and (if it was already listed for infection) before the
				// previously listed infection event,
                // and if the infection time is at the latest at t* + T, then list it
				if ((t <= n[me].time) && (t < n[you].time) && (t <= (g.true_start + g.T))) {
					// set you's infection time
                    n[you].time = t;
                    // store who you gets infected by
                    n[you].inf_by = me;
                    // if not listed before, then extend the heap
					if (n[you].heap == NONE) {
						g.heap[++g.nheap] = you;
						n[you].heap = g.nheap;
					}
                    // this works bcoz the only heap relationship that can be violated is the one between you and its parent
					up_heap(n[you].heap);
				}
			}
		}
	}

    // add me to g.s
	g.s[g.ns++] = me;
    
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// this routine runs one SIR outbreak from a random starting node

void sir (unsigned int source, unsigned int start) {
    
    // declare integers i
	unsigned int i;
	
    // initialize size of outbreak to 0
	g.ns = 0;
	
	// infect the source at time start (0, ..., g.dur)
	n[source].time = start;
	
    // set source on position 1 of heap
    n[source].heap = 1;
    
    // add source to heap
	g.heap[g.nheap = 1] = source;

	// run the outbreak
	while (g.nheap) infect();

	// clean
	for (i = 0; i < g.ns; i++) n[g.s[i]].heap = n[g.s[i]].time = NONE;
    
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// function to generate ground truth

void ground_truth () {
    
    // declare unsigned integer
	unsigned int i, upper, lower;
    
    // upper and lower limit of starting times
    upper = g.dur - g.T; lower = (NUMT - 1) * g.T + 1;
    
    // set all true_state values to 0 (= susceptible)
    // true_state describes the disease status at g.dur
    for (i = 0; i < g.n; i++) n[i].true_state = 0;
    
    // simulate until we get outbreak that has at least size of MINOUTBREAK
    do {
        
        // sample starting time
        g.true_start = (pcg_32() % (upper - lower + 1)) + lower;
                
        // sample true source index
        g.true_source_idx = pcg_32_bounded(g.nps);
        
        // store true source
        g.true_source = g.ps[g.true_source_idx];
        
        // run sir() process
        sir(g.true_source, g.true_start);
        
        // if outbreak size is smaller than MINOUTBREAK, set inf_time and rec_time back to NONE
        if (g.ns < MINOUTBREAK) for (i = 0; i < g.ns; i++) n[g.s[i]].inf_time = n[g.s[i]].rec_time = n[g.s[i]].inf_by = NONE;
        
	} while (g.ns < MINOUTBREAK); 
    
    // print info to console
    printf("True seed node is %d (getting infected at time %d)\n", g.true_source, g.true_start);
    printf("Number of infected or recovered nodes: %d\n", g.ns);
    printf("First observation of outbreak occurrs at %d\n", g.true_start + g.T);
    
    // store ground truth   
    // loop over all activated nodes in g.s
    for (i = 0; i < g.ns; i++) {
        // if the recovery time is later than t* + T, the state is 1 (= infectious)
        if (n[g.s[i]].rec_time > (g.true_start + g.T)) {
            n[g.s[i]].true_state = 1;
        } else {
            // otherwise, it is 2 (= recovered)
            n[g.s[i]].true_state = 2;
        }
        printf("Node %d: Infected by %d at %d. Recovery time is %d, hence the true state is %d\n", g.s[i], n[g.s[i]].inf_by, n[g.s[i]].inf_time, n[g.s[i]].rec_time, n[g.s[i]].true_state);
    }
    
    // store size of ground truth and time of inference in GLOBALS
    g.ns_gt = g.ns;
    g.t_now = g.true_start + g.T;
    
    // set inf_time and rec_time of nodes back to NONE
    for (i = 0; i < g.ns; i++) n[g.s[i]].inf_time = n[g.s[i]].rec_time = n[g.s[i]].inf_by = NONE;
    
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// function to create simulation results for all possible start configurations

void simulate (unsigned int nexp, INF *inf, PRIOR *pr, OUTF *outf) {
    
    // declare unsigned integers
    unsigned int i, j, k, T, start;
    
    // declare start and end for storing CPU run times
    clock_t t0, t1;
    
    // record start time t0
    t0 = clock();
    
    // compute earliest possible start
    start = g.t_now - pr->n_durs;
    
    // find all possible sources with deterministic SIR
    find_active_sources(inf, start);
    
    // allocate memory to array with MCS probabilities
    p = calloc(inf->n_sources, sizeof(MCS));
    
    // ------------------------------------------------
    // run Monte-Carlo simulations and store results
    
    // 1. Loop over sources and memory allocation
    for (i = 0; i < inf->n_sources; i++) {
        
        // allocate memory to outer arrays
        p[i].inf = malloc(g.n * sizeof(double *));
        p[i].rec = malloc(g.n * sizeof(double *));
        
        // loop over nodes in graph g.n
        for (j = 0; j < g.n; j++) {
            
            // allocate memory to inner arrays
            // for each node we want to store the prob. of 
            // being infected / recovered if the current source 
            // started the infection at possible starting times
            p[i].inf[j] = calloc(pr->n_durs, sizeof(double));
            p[i].rec[j] = calloc(pr->n_durs, sizeof(double));
            
        }
        
        // 2. Loop over possible values for T
        for (T = 0; T < pr->n_durs; T++) {
            
            // compute epidemic start
            start = g.t_now - T;
            
            // 3. Simulate NSIM times
            for (j = 0; j < NSIM; j++) {
                
                // run SIR from source at index i, starting at time t
                sir(inf->sources[i], start);
                
                // loop over all activated nodes in g.s
                for (k = 0; k < g.ns; k++) {
                    
                    // store results for each activated node
                    if (n[g.s[k]].rec_time > g.t_now) {
                        // add 1.0 to inf. array if node's recovery time is later than g.dur
                        p[i].inf[g.s[k]][T] += 1.0;
                    } else {
                        // add 1.0 to rec. array if node's recovery time is before g.dur
                        p[i].rec[g.s[k]][T] += 1.0;
                    }
            
                    // set inf_time and rec_time of nodes back to NONE
                    n[g.s[k]].inf_time = n[g.s[k]].rec_time = n[g.s[k]].inf_by = NONE;
                    
                }
                
            }
            
            // 4. Loop over nodes in graph and normalize in order to get probabilities
            for (j = 0; j < g.n; j++) {
            
                // if element in array is larger than 0.0, normalize by NSIM
                if (p[i].inf[j][T] > 0.0) p[i].inf[j][T] /= NSIM;
                if (p[i].rec[j][T] > 0.0) p[i].rec[j][T] /= NSIM;
            
            }
            
        }
  
    }
    
    // record end time t1
    t1 = clock();
    
    // store simulation run time
    outf->sim_time[nexp] = ((double)(t1 - t0)) / CLOCKS_PER_SEC;
 
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -