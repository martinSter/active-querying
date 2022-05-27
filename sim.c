// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// simulation functions

// include header file
#include "run.h"

// declare external variables
extern GLOBALS g;
extern NODE *n;
extern MCS *p;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// this routine does the bookkeeping for an infection event

void infect () {
    
    // declare unsigned integers
	unsigned int i, you, me = g.heap[1];
    
    // declare floats
	float t, now = n[me].time;

    // take the newly infected off the heap
	del_root();
    
    // store infection time of me
    n[me].inf_time = n[me].time;
    
	// get the recovery time of me
    // multiply by g.beta to get the recovery time (which is 1 by default)
	n[me].time += g.rexp[pcg_16()] * g.beta;
    
    // store recovery time of me
    n[me].rec_time = n[me].time;

	// go through the neighbors of the infected node
	for (i = 0; i < n[me].deg; i++) {
        
        // get neighbor at index i
		you = n[me].nb[i];
        
        // if you is S, you can be infected
		if (S(you)) {
        
            // get the infection time of you
			t = now + g.rexp[pcg_16()];
            
            // make sure that t is before recovery time of me and...
            // t is before current infection time of you and...
            // t is smaller or equal to TEND
			if ((t < n[me].time) && (t < n[you].time) && (t <= TEND)) {
                
                // if that is satisfied, set time of you to t
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
    
    // add me to g.s
	g.s[g.ns++] = me;
    
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// this routine runs one SIR outbreak from a random starting node

void sir (unsigned int source) {
    
    // declare integers i
	unsigned int i;
	
    // initialize size of outbreak to 0
	g.ns = 0;
	
	// infect the source at time 0
	n[source].time = 0.0;
	
    // set source on position 1 of heap
    n[source].heap = 1;
    
    // add source to heap
	g.heap[g.nheap = 1] = source;

	// run the outbreak
	while (g.nheap) infect();

	// clean
	for (i = 0; i < g.ns; i++) {
        n[g.s[i]].heap = NONE;
        n[g.s[i]].time = DBL_MAX;
    }
    
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// function to generate ground truth

void ground_truth () {
    
    // declare unsigned integers
	unsigned int i, j;
    
    // declare floats
    float min, max, now;
    
    // allocate memory to ground_truth
    for (i = 0; i < g.n; i++) n[i].ground_truth = calloc(g.ntinf, sizeof(unsigned int));
    
    // simulate until we get outbreak that has at least size of MINOUTBREAK
    do {
        
        // store true source
        g.true_source = pcg_32_bounded(g.n);
        
        // run sir() process
        sir(g.true_source);
        
        // if outbreak size is smaller than MINOUTBREAK, set inf_time and rec_time back to DBL_MAX
        if (g.ns < MINOUTBREAK) {
            for (i = 0; i < g.ns; i++) {
                n[g.s[i]].inf_time = n[g.s[i]].rec_time = DBL_MAX;
                n[g.s[i]].inf_by = NONE;
            }
        }
        
	} while (g.ns < MINOUTBREAK);
    
    // print info to console
    printf("True seed node is %d\n", g.true_source);
    printf("Number of infected or recovered nodes: %d\n", g.ns);
    
    // store ground truth   
    // loop over all activated nodes in g.s
    for (i = 0; i < g.ns; i++) {
        // loop over inference times
        for (j = 0; j < g.ntinf; j++) {
            // set min to infection time of node
            min = n[g.s[i]].inf_time;
            // set max to recovery time of node
            max = n[g.s[i]].rec_time;
            // set now to current inference time
            now = g.tinf[j];
            // if now is within infectious period of node, set state to 1 (infectious)
            if ((now - min) * (now - max) < 0) n[g.s[i]].ground_truth[j] = 1;
            // if now is larger or equal to max, set state to 2 (recovered)
            if (now >= max) n[g.s[i]].ground_truth[j] = 2;
        }
    }
    
    // store size of ground truth in GLOBALS
    g.ns_gt = g.ns;
    
    // set inf_time and rec_time of nodes back to DBL_MAX
    for (i = 0; i < g.ns; i++) {
        n[g.s[i]].inf_time = n[g.s[i]].rec_time = DBL_MAX;
        n[g.s[i]].inf_by = NONE;
    }
    
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// function to create simulation results for all possible sources

void monte_carlo () {
    
    // declare unsigned integers
    unsigned int i, j, k, h, os, r0;
    
    // declare integers
    float t, min, max, now;
    
    // initialize g.os and g.r0
    g.os = 0.0; g.r0 = 0.0;
    
    // ------------------------------------------------
    // store inference times
    
    // set number of inference times
    g.ntinf = (TEND - 0.0) / TSTEPS;
    
    // allocate memory
    g.tinf = malloc(g.ntinf * sizeof(unsigned int));
    
    // set t to 0 (start of time period)
    t = 0.0;
    
    // loop over inference times
    for (i = 0; i < g.ntinf; i++) {
        // compute current inference time
        t += TSTEPS;
        // store inference time in g.tinf
        g.tinf[i] = t;
    }
    
    // ------------------------------------------------
    // allocate memory to p and all its subarrays
    
    // allocate memory to array with MCS probabilities
    p = calloc(g.n, sizeof(MCS));
    
    // loop over sources
    for (i = 0; i < g.n; i++) {
        // allocate memory to outer arrays
        p[i].inf = malloc(g.n * sizeof(double *));
        p[i].rec = malloc(g.n * sizeof(double *));
        // loop over nodes in graph g.n
        for (j = 0; j < g.n; j++) {
            // allocate memory to inner arrays
            p[i].inf[j] = calloc(g.ntinf, sizeof(double));
            p[i].rec[j] = calloc(g.ntinf, sizeof(double));
        }
    }
    
    // ------------------------------------------------
    // run Monte-Carlo simulations and store results
    
    // loop over sources
    for (i = 0; i < g.n; i++) {
        
        // initialize os to 0
        os = 0; r0 = 0;
        
        // simulate NSIM times
        for (j = 0; j < NSIM; j++) {
            
            // run SIR from source i
            sir(i);
            
            // store outbreak size
            os += g.ns;
            
            // loop over all activated nodes in g.s
            for (k = 0; k < g.ns; k++) {
                
                // if node was infected by source, then increment r0 by 1
                if (n[g.s[k]].inf_by == i) r0++;
                
                // loop over inference times
                for (h = 0; h < g.ntinf; h++) {
                    
                    // set min to infection time of node
                    min = n[g.s[k]].inf_time;
                    // set max to recovery time of node
                    max = n[g.s[k]].rec_time;
                    // set now to current inference time
                    now = g.tinf[h];
                    
                    // if now is within infectious period of node, add 1.0 to infectious array
                    if ((now - min) * (now - max) < 0) p[i].inf[g.s[k]][h] += 1.0;
                    
                    // if now is larger or equal to max, add 1.0 to recovered array
                    if (now >= max) p[i].rec[g.s[k]][h] += 1.0;
                    
                }
        
                // set inf_time and rec_time of nodes back to DBL_MAX
                n[g.s[k]].inf_time = n[g.s[k]].rec_time = DBL_MAX;
                
                // set inf_by back to NONE
                n[g.s[k]].inf_by = NONE;
                
            }
            
        }
        
        // store average outbreak size of source i in g.os
        g.os += ((double) os / NSIM);
        g.r0 += ((double) r0 / NSIM);
        
    }
    
    // print average outbreak size and R0 to console
    printf("Average outbreak size: %f\n", g.os / g.n);
    printf("Average R0: %f\n", g.r0 / g.n);
    
    // ------------------------------------------------
    // normalize in order to get probabilities
    
    // loop over sources
    for (i = 0; i < g.n; i++) {
        // loop over nodes in graph
        for (j = 0; j < g.n; j++) {
            // loop over inference times
            for (k = 0; k < g.ntinf; k++) {
                // add-one smoothing and normalization
                p[i].inf[j][k] = (p[i].inf[j][k] + 1) / (NSIM + 1);
                p[i].rec[j][k] = (p[i].rec[j][k] + 1) / (NSIM + 1);
                // if element in array is larger than 0.0, normalize by NSIM
                //if (p[i].inf[j][k] > 0.0) p[i].inf[j][k] /= (NSIM + 1);
                //if (p[i].rec[j][k] > 0.0) p[i].rec[j][k] /= (NSIM + 1);
            }
            
        }
        
    }
    
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -