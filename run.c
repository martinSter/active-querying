// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Run the active learning approach for epidemics on temporal networks
// - Note: the simulation part of the approach is based on code written by 
//         Petter Holme (https://github.com/pholme/tsir)

// include header file
#include "run.h"

// declare g as a GLOBALS struct
GLOBALS g;

// declare array of NODE structs
NODE *n;

// declare array with Monte-Carlo probs
MCS *p;

// declare inference and prior struct
INF inf; PRIOR pr;

// declare output structs
OUTF outf;
OUTT outt_maxp; OUTT outt_ucty; OUTT outt_akld; OUTT outt_onehop; OUTT outt_randq;
OUTQ out_maxp; OUTQ out_ucty; OUTQ out_akld; OUTQ out_onehop; OUTQ out_randq;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// main function

int main (int argc, char *argv[]) {
    
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    // Some declarations
    
    // declare unsigned integers
	unsigned int i, j, k;
    
    // declare some doubles
	double d, x;
    
    // declare pointer to FILE
	FILE *fp;
	
	// just a help message
	if ((argc < 7) || (argc > 7)) {
		fprintf(stderr, "usage: ./run [data/file] [beta] [mu] [T] [undirected = 0, directed = 1] [seed]\n");
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
    
    // call read_data function
	read_data(fp, atoi(argv[5]));
    
    // close file
	fclose(fp);
    
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    // Setup
   
    // add infection probability beta to GLOBALS
    g.beta = atof(argv[2]);
	
    // add recovery probability beta to GLOBALS
    g.mu = atof(argv[3]);
    
    // add T (time until first obs.) to GLOBALS
    g.T = atoi(argv[4]);
    
    // initialize random number generator
    g.state = (uint64_t) strtoull(argv[6], NULL, 10);
    
    // compute 1 / log(1 - beta) --> inverse probability transform
	d = 1.0 / log(1.0 - g.beta);
	
    // create a large number (65536) of random numbers created according 
    // to geometric distribution (discrete time!)
    // note that 2^16 = 65536
    for (i = 0; i < 0x10000; i++) {
        // inverse probability transform
        x = d * log((i + 1) / 65536.0);
		// add random samples to array
        g.rnd2inx[i] = (x > USHRT_MAX) ? USHRT_MAX : x;
	}
    
    // run the random number creator a first time
    // (otherwise it outputs 0 the first time it is used... weird)
    pcg_32_bounded(0x10000);

	// allocating the heap (N + 1) because its indices are 1,...,N
    // g.heap is a pointer, so here we assign it the address where it will point to
	g.heap = malloc((g.n + 1) * sizeof(unsigned int));
    
    // allocate memory to g.s (array containing all nodes that get infected)
    g.s = calloc(g.n, sizeof(unsigned int));
    
    // initialize
	for (i = 0; i < g.n; i++) {
        n[i].heap = n[i].time = NONE;
        n[i].inf_time = n[i].rec_time = n[i].inf_by = NONE;
    }
    
    // find all possible sources (nodes with out-edges)
    // we will sample from this to generate ground truth
    find_possible_sources();
    
    // initialize the prior over starting times
    // type = 0 (geometric), type = 1 (delta)
    init_prior(&pr, 0);
    
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    // Allocate memory to output structs
    
    // allocate memory to arrays in outf
    outf.nactives = malloc(NEXP * sizeof(unsigned int));
    outf.nsources = malloc(NEXP * sizeof(unsigned int));
    outf.true_sources = malloc(NEXP * sizeof(unsigned int));
    outf.ranks = malloc(NEXP * sizeof(unsigned int));
    outf.posteriors = malloc(NEXP * sizeof(double));
    outf.css = malloc(NEXP * sizeof(unsigned int));
    outf.sim_time = malloc(NEXP * sizeof(double));
    
    // allocate memory to arrays in outt
    outt_maxp.run_time = malloc(NEXP * sizeof(double));
    outt_ucty.run_time = malloc(NEXP * sizeof(double));
    outt_akld.run_time = malloc(NEXP * sizeof(double));
    outt_onehop.run_time = malloc(NEXP * sizeof(double));
    outt_randq.run_time = malloc(NEXP * sizeof(double));
    
    // set up output structs for querying
    set_up(&out_maxp);
    set_up(&out_ucty);
    set_up(&out_akld);
    set_up(&out_onehop);
    set_up(&out_randq);
    
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    // Experiments
    
    // loop over experiments
    for (i = 0; i < NEXP; i++) {
        
        // print line to console
        printf("--------------------------------------------\n");
        
        // print number of experiment to console
        printf("NUMBER OF EXPERIMENT: %d\n", i + 1);
        
        // generate ground truth outbreak
        ground_truth();
        
        // draw first observed node (for querying part of the experiment)
        // need to do that here because we need this to determine the 
        // sources that are possible (every node that can reach first obs. node)
        draw_first_node();
        
        // print line to console
        printf("--------------------------------------------\n");
        
        // Monte Carlo simulations
        simulate(i, &inf, &pr, &outf);
        
        // run FULL inference
        full_inference(i, &inf, &pr, &outf);
        
        // print line to console
        printf("--------------------------------------------\n");
        
        // run QUERYING
        infer_query(i, &inf, &pr, maxp, &out_maxp, &outt_maxp, 0);
        infer_query(i, &inf, &pr, ucty, &out_ucty, &outt_ucty, 0);
        infer_query(i, &inf, &pr, akld, &out_akld, &outt_akld, 0);
        infer_query(i, &inf, &pr, onehop, &out_onehop, &outt_onehop, 0);
        infer_query(i, &inf, &pr, randq, &out_randq, &outt_randq, 0);
        
        // cleaning up p
        for (k = 0; k < inf.n_sources; k++) {
            for (j = 0; j < g.n; j++) {
                // free inner arrays
                free(p[k].inf[j]);
                free(p[k].rec[j]);
            }    
            // free outer arrays
            free(p[k].inf);
            free(p[k].rec);
        }
        
        // info to console
        printf("Finished cleaning up...\n");
    
    }
    
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    // Write results to file
    
    // export results from FULL inference experiments
    export_output(&outf);
    
    // export results from QUERYING experiments
    export_output_querying(&out_maxp, "maxp");
    export_output_querying(&out_ucty, "ucty");
    export_output_querying(&out_akld, "akld");
    export_output_querying(&out_onehop, "onehop");
    export_output_querying(&out_randq, "randq");
    
    // export run times
    export_run_times(&outt_maxp, &outt_ucty, &outt_akld, &outt_onehop, &outt_randq);
    
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    // Cleaning up
    
    // print info to console
    printf("Cleaning up...\n");
        
	// cleaning up
	for (i = 0; i < g.n; i++) {
		for (j = 0; j < n[i].deg; j++) free(n[i].t[j]);
        // free all arrays in NODE struct n
		free(n[i].nb);
		free(n[i].nc);
		free(n[i].t);
	}
    
    // free array n of NODE structs and heap and ground_truth
	free(n); free(g.heap); free(g.ground_truth); free(g.s); free(g.ps); free(g.tinf);
    
    // free outf arrays
    free(outf.css); free(outf.posteriors); free(outf.ranks); 
    free(outf.true_sources); free(outf.nsources); free(outf.nactives); free(outf.sim_time);
    
    // free outt arrays
    free(outt_maxp.run_time); free(outt_ucty.run_time); free(outt_akld.run_time);
    free(outt_onehop.run_time); free(outt_randq.run_time);
    
    // clean up output structs from querying
    clean_up(&out_maxp);
    clean_up(&out_ucty);
    clean_up(&out_akld);
    clean_up(&out_onehop);
    clean_up(&out_randq);
    
    // free memory allocated to sources and prior structs
    free(inf.sources); free(pr.priors);
    
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	return 0;
    
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
