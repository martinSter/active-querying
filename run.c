// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Run the active learning approach for epidemics on temporal networks
// - Note: the simulation part of the approach is based on beautiful code 
//         written by Petter Holme (https://github.com/pholme/tsir)

// include header file
#include "run.h"

// declare g as a GLOBALS struct
GLOBALS g;

// declare array of NODE structs
NODE *n;

// declare array with Monte-Carlo probs
MCS *p;

// declare struct for posterior distributions
POSTERIOR post;

// declare output structs
OUTPUT out_full;
OUTPUTQ out_maxp;
OUTPUTQ out_ucty;
OUTPUTQ out_akld;
OUTPUTQ out_onehop;
OUTPUTQ out_randq;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Main function

int main (int argc, char *argv[]) {
    
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    // Some declarations and run command control
    
    // declare unsigned integers
	unsigned int i, j;
    
    // declare pointer to FILE
	FILE *fp;
	
	// just a help message
    // note that mu = 1 (recovery rate) is default
	if ((argc < 4) || (argc > 4)) {
		fprintf(stderr, "usage: ./run [data/file] [beta] [seed]\n");
		return 1;
	}
    
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    // Import the network from file

	// open network data file
	fp = fopen(argv[1], "r");
	if (!fp) {
        // print error if we cannot open file
		fprintf(stderr, "can't open '%s'\n", argv[1]);
		return 1;
	}
    
    // call read_data function to import network data
	read_data(fp);
    
    // close file
	fclose(fp);
    
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    // Setup
    
    // add infection probability beta to GLOBALS
    g.beta = atof(argv[2]);
    
    // for a quick check, output transmission rate
    printf("Transmission rate: %3.2f\n", g.beta);
    
    // initialize random number generator
    g.state = (uint64_t) strtoull(argv[3], NULL, 10);
	
    // create a large number (65536) of random numbers from exponential distr. (time until infection)   
    for (i = 0; i < 0x10000; i++)
		g.rexp[i] = -log((i + 1.0) / 0x10000) / g.beta;

	// allocating the heap (N + 1) because its indices are 1,...,N
    // g.heap is a pointer, so here we assign it the address where it will point to
	g.heap = malloc((g.n + 1) * sizeof(unsigned int));
    
    // allocate memory to g.s (array containing all nodes that get infected)
    g.s = calloc(g.n, sizeof(unsigned int));
    
    // initialize
	for (i = 0; i < g.n; i++) {
        // NONE and DBL_MAX are defined in run.h
        n[i].heap = n[i].inf_by = NONE;
        n[i].time = n[i].inf_time = n[i].rec_time = DBL_MAX;
    }
    
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    // Monte-Carlo simulations from all possible sources
    
    // run Monte-Carlo simulations from all nodes
    // monte_carlo() is defined in sim.c
    monte_carlo();
    
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    // Set up output structs
    
    // allocate memory to arrays in out_full
    out_full.truesource = malloc(NEXP * sizeof(unsigned int));
    out_full.qmap = malloc(NEXP * sizeof(unsigned int *));
    out_full.qrand = malloc(NEXP * sizeof(unsigned int *));
    out_full.ranks = malloc(NEXP * sizeof(unsigned int *));
    out_full.nsources = malloc(NEXP * sizeof(unsigned int *));
    out_full.nactive = malloc(NEXP * sizeof(unsigned int *));
    
    // allocate memory to inner arrays
    for (i = 0; i < NEXP; i++) {
        out_full.qmap[i] = malloc(g.ntinf * sizeof(unsigned int));
        out_full.qrand[i] = malloc(g.ntinf * sizeof(unsigned int));
        out_full.ranks[i] = malloc(g.ntinf * sizeof(unsigned int));
        out_full.nsources[i] = malloc(g.ntinf * sizeof(unsigned int));
        out_full.nactive[i] = malloc(g.ntinf * sizeof(unsigned int));
    }
    
    // set up output structs for querying
    // set_up() is defined in misc.c
    set_up(&out_maxp);
    set_up(&out_ucty);
    set_up(&out_akld);
    set_up(&out_onehop);
    set_up(&out_randq);
    
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    // Experiments
    
    // loop from 0 to NEXP
    for (i = 0; i < NEXP; i++) {
        
        // print line to console
        printf("--------------------------------------------\n");
        
        // print number of experiment to console
        printf("NUMBER OF EXPERIMENT: %d\n", i + 1);
        
        // generate ground truth outbreak
        // ground_truth() is defined in sim.c
        ground_truth();
        
        // add true source to OUTPUT out_full struct
        out_full.truesource[i] = g.true_source;
        
        // print line to console
        printf("--------------------------------------------\n");
        
        // run FULL inference
        // full_inference() is defined in inf.c
        full_inference(i, &post, &out_full);
        
        // run QUERYING
        // run_querying() is defined in inf.c
        run_querying(i, &post, &out_maxp, &out_ucty, &out_akld, &out_onehop, &out_randq);
        
        // clean up ground_truth
        for (j = 0; j < g.n; j++) free(n[j].ground_truth);
    
    }
    
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    // Write results to file
    
    // export results from FULL inference experiments
    export_output(&out_full);
    
    // export results from QUERYING experiments
    export_output_querying(&out_maxp, "maxp");
    export_output_querying(&out_ucty, "ucty");
    export_output_querying(&out_akld, "akld");
    export_output_querying(&out_onehop, "onehop");
    export_output_querying(&out_randq, "randq");
    
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    // Cleaning up
    
    // print info to console
    printf("Cleaning up...\n");
        
	// cleaning up
	for (i = 0; i < g.n; i++) free(n[i].nb);
    
    // free array n of NODE structs and heap and ground_truth
	free(n); free(g.heap); free(g.ground_truth); free(g.s); free(g.tinf);
    
    // free node probs.
    for (i = 0; i < g.n; i++) {
		for (j = 0; j < g.n; j++) {
            // free inner arrays
            free(p[i].inf[j]);
            free(p[i].rec[j]);
        }    
        // free outer arrays
		free(p[i].inf);
        free(p[i].rec);
    }
    
    // free ranks in out_full
    for (i = 0; i < NEXP; i++) {
        free(out_full.qmap[i]);
        free(out_full.qrand[i]);
        free(out_full.ranks[i]);
        free(out_full.nsources[i]);
        free(out_full.nactive[i]);
    }
    
    // free full arrays
    free(out_full.truesource); free(out_full.qmap); free(out_full.qrand); free(out_full.ranks); free(out_full.nsources); free(out_full.nactive);
    
    // clean up output structs from querying
    // clean_up() is defined in misc.c
    clean_up(&out_maxp);
    clean_up(&out_ucty);
    clean_up(&out_akld);
    clean_up(&out_onehop);
    clean_up(&out_randq);
    
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	return 0;
    
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
