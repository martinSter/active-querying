// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// miscellaneous routines

// part of it is code for temporal network SIR by Petter Holme (2018)

// include header file
#include "run.h"

// declare external variables
extern GLOBALS g;
extern NODE *n;
extern MCS *p;

// declare external variables
extern EVID evid_random;
extern EVID evid_onehop;
extern EVID evid_maxp;
extern EVID evid_hard;
extern EVID evid_soft;
extern EVID evid_esr;
extern EVID evid_kl;

// declare and define pointer 'alloc'
unsigned int *alloc;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// reads the network, assumes an edge list with vertex label 0,N-1
// if your network has nodes with degree zero, make sure that none of them is
// the node with largest index

void read_data (FILE *fp) {
    
    // declare unsigned integers
	unsigned int i, me, you;
    
    // initialize number of nodes
	g.n = 0; g.e = 0;

	// scan the system size to find g.n
	while (2 == fscanf(fp, "%u %u\n", &me, &you)) {
		if (g.n < me) g.n = me;
		if (g.n < you) g.n = you;
        // increment number of edges by 1
        g.e++;
	}

    // add 1 to g.n since first node is 0
	g.n++;
    
    // prints out number of nodes and number of edges
    printf("Number of nodes: %d\nNumber of edges: %d\n", g.n, g.e);

    // allocate memory to n
	n = calloc(g.n, sizeof(NODE));

    // rewind file
	rewind(fp);

	// scan the degrees
	while (2 == fscanf(fp, "%u %u\n", &me, &you)) {
		n[me].deg++;
		n[you].deg++;
	}

	// allocate adjacency lists
	for (i = 0; i < g.n; i++) {
        // allocate memory to nb arrays
		n[i].nb = malloc(n[i].deg * sizeof(unsigned int));
        // reset deg to 0 as we need it to allocate neighbors
		n[i].deg = 0;
	}

    // rewind file
	rewind(fp);

	// fill adjacency lists
	while (2 == fscanf(fp, "%u %u\n", &me, &you)) {
		n[me].nb[n[me].deg++] = you;
		n[you].nb[n[you].deg++] = me;
	}
    
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// set up output structs for querying

void set_up (OUTPUTQ *out) {
    
    // declare integers
    unsigned int i, j;
    
    // allocate memory to arrays in out (first level)
    out->qmap = malloc(NEXP * sizeof(unsigned int **));
    out->ranks = malloc(NEXP * sizeof(unsigned int **));
    out->nsources = malloc(NEXP * sizeof(unsigned int **));
    
    // loop over number of experiments
    for (i = 0; i < NEXP; i++) {
        
        // allocate memory to arrays on second level
        out->qmap[i] = malloc(g.ntinf * sizeof(unsigned int *));
        out->ranks[i] = malloc(g.ntinf * sizeof(unsigned int *));
        out->nsources[i] = malloc(g.ntinf * sizeof(unsigned int *));
        
        // loop over inference times
        for (j = 0; j < g.ntinf; j++) {
            
            // allocate memory to array on third level
            out->qmap[i][j] = malloc(NQUERIES * sizeof(unsigned int));
            out->ranks[i][j] = malloc(NQUERIES * sizeof(unsigned int));
            out->nsources[i][j] = malloc(NQUERIES * sizeof(unsigned int));
            
        }
        
    }
    
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// clean up output structs for querying

void clean_up (OUTPUTQ *out) {
    
    // declare integers
    unsigned int i, j;
    
    // loop over number of experiments
    for (i = 0; i < NEXP; i++) {
        
        // loop over inference times
        for (j = 0; j < g.ntinf; j++) {
            
            // free arrays on third level
            free(out->qmap[i][j]);
            free(out->ranks[i][j]);
            free(out->nsources[i][j]);
            
        }
            
        // free arrays on second level
        free(out->qmap[i]);
        free(out->ranks[i]);
        free(out->nsources[i]);
        
    }
    
    // free arrays on first level
    free(out->qmap); free(out->ranks); free(out->nsources);
    
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// drawing the first observed infected node

void draw_first_node (unsigned int tinf_idx) {
    
    // declare integers
    unsigned int i, active = 0;
    
    // create array for active nodes
    unsigned int active_nodes[g.n];
    
    // loop over all nodes in graph
    for (i = 0; i < g.n; i++) {
        // check if node is activated
        if (n[i].ground_truth[tinf_idx] > 0) {
            // add node i to active_nodes
            active_nodes[active] = i;
            // increment number of active nodes by 1
            active++;
        }
    }
    
    // random draw of integer in active nodes
    i = pcg_32_bounded(active);
    
    // store first observed
    g.fin = active_nodes[i];
    
    // print to console
    // printf("First observed node is %u\n", g.fin);
    
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// initialize evidence with first observed node

void initialize_evidence (EVID *evid) {
    
    // declare integer
    unsigned int i;
    
    // set e to 0 for each node
    for (i = 0; i < g.n; i++) n[i].e = 0;
       
    // set number of evidence nodes to 0
    evid->n_evid = 0;
    
    // allocate memory to array of nodes in evidence
    evid->nodes = calloc(1, sizeof(unsigned int));
    
    // add first observed node to evidence
    evid->nodes[evid->n_evid++] = g.fin;
    
    // set e to 1 for g.fin
    n[g.fin].e = 1;
    
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// update evidence

void update_evidence (unsigned int v, EVID *evid) {
    
    // reallocate memory
    evid->nodes = realloc(evid->nodes, (evid->n_evid + 1) * sizeof(unsigned int));
    
    // add new node v to evidence
    evid->nodes[evid->n_evid++] = v;
    
    // set e to 1 for new node
    n[v].e = 1;

}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// find max. node in array of doubles

unsigned int find_max (double array[]) {
    
    // declare integers
    unsigned int i, nmax = 0;
    
    // declare double
    double max = 0.0;
    
    // array for max. nodes
    unsigned int pn[g.n];
    
    // find node where element is maximal
    for (i = 0; i < g.n; i++) {
        
        // jump to next iteration if array element is 0
        if (array[i] == 0.0) continue;
        
        // check if current array element is larger than current max
        if (array[i] > max) {
            
            // set max to current array element
            max = array[i];
            
            // set nmax back to 0
            nmax = 0;
        }
        
        // add i-th element to pn if it corresponds to max. value
        if(array[i] == max) pn[nmax++] = i;

    }
        
    // return node
    if (nmax == 1) return pn[0];
    else if (nmax > 1) return pn[pcg_32_bounded(nmax)];
    else return NONE;  

}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// sort log-likelihoods and find rank of true source

void sort_rank (POSTERIOR *post) {
    
    // declare integers
    unsigned int i, j, temp1, rank, nn, ns;
    
    // declare doubles
    double min, temp2, rank_final = -1.0;
    
    // initialize arrays on the stack
    unsigned int source_indices[g.n]; // this is needed for the sorting
    double log_liks[g.n]; // this is needed for the sorting
    double ranks[g.n]; // ranks can be decimal numbers
    
    // *******************************************************************
    // Find min. log-likelihood, then set all invalid sources to min - 1.0
    
    // initialize min to first element in log_liks
    min = post->log_liks[0];
    
    // initialize ns to number of sources
    ns = g.n;
    
    // loop over all sources
    for (i = 0; i < g.n; i++) {
        
        // find min. element in log_liks
        if (post->log_liks[i] < min) min = post->log_liks[i];
        
        // reduce ns by 1 if source is invalid
        if (post->log_liks[i] == INVALID) ns--;
     
    }
    
    // ****************************************
    // Sort log-likelihoods in descending order

    // loop over all sources
    for (i = 0; i < g.n; i++) {
        // fill source_indices with i
        source_indices[i] = i;
        // fill log_liks
        log_liks[i] = (post->log_liks[i] == INVALID) ? (min - 1.0) : post->log_liks[i];
    }
 
    // loop over array elements
    for (i = 0; i < g.n; i++) {
        
        // loop over all array elements after i-th element
        for (j = i + 1; j < g.n; j++) {
            
            // swap elements if i-th element is smaller than j-th element
            if (log_liks[i] < log_liks[j]) {
                
                // store source and lik at i in temp
                temp1 = source_indices[i];
                temp2 = log_liks[i];
                
                // set i to j
                source_indices[i] = source_indices[j];
                log_liks[i] = log_liks[j];
                
                // set j from temp
                source_indices[j] = temp1;
                log_liks[j] = temp2;
            }
        }
    }
    
    // store max. log-likelihood in GLOBALS
    g.maxlog = log_liks[0];
    
    // store qmap
    g.qmap = source_indices[0];
    
    // ******************************************
    // Compute ranks and find rank of true source
    
    // initialize values
    i = 0; rank = 1;
    
    // loop over array
    while (i < g.n) {
        // set j to current value of i
        j = i;
        // count how many elements with same log-liks
        while (j < g.n - 1 && log_liks[j] == log_liks[j + 1]) j++;
        // set nn
        nn = j - i + 1;
        // set ranks
        for (j = 0; j < nn; j++) ranks[j + i] = rank + (nn - 1) * 0.5;
        // increment rank and i
        rank += nn;
        i += nn;
    }
    
    // find rank of true source
    for (i = 0; i < g.n; i++) if (source_indices[i] == g.true_source) rank_final = ranks[i];
    
    // store results in GLOBALS
    g.rank = rank_final;
    g.nsources = ns;
    
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// randomly select node from all activated nodes

unsigned int random_draw (unsigned int inf_time_idx) {
    
    // declare integers
    unsigned int rn, i, j = 0;
    
    // create array
    unsigned int possible_nodes[g.n];
    
    // add all activated nodes to array
    for (i = 0; i < g.n; i++) if (n[i].ground_truth[inf_time_idx] > 0) possible_nodes[j++] = i;
        
    // randomly pick one
    rn = possible_nodes[pcg_32_bounded(j)];
        
    // return node
    return rn;
    
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// write results to file

void export_output (OUTPUT *out) {
    
    // declare unsigned integers
    unsigned int i, j;
    
    // declare pointer to FILE
	FILE *fp;
    
    // open file
	fp = fopen("output/output_FULL_true_sources.txt", "w");
    if (!fp) {
        // print error if we cannot open file
		fprintf(stderr, "can't open file\n");
		return;
	}
    
    // loop over experiments and write results to file
    for (i = 0; i < NEXP; i++) fprintf(fp, "%u\n", out->truesource[i]);

    // close file
	fclose(fp);
    
    // --------------------------------------------------
    
    // open file
	fp = fopen("output/output_FULL_qmap.txt", "w");
    if (!fp) {
        // print error if we cannot open file
		fprintf(stderr, "can't open file\n");
		return;
	}
    
    // loop over experiments
    for (i = 0; i < NEXP; i++) {
        
        // loop over inference times
        for (j = 0; j < g.ntinf; j++) {
            
            // write results to file
            fprintf(fp, "%u", out->qmap[i][j]);
            
            // either ";" or newline after a value
            if (j == (g.ntinf - 1)) fprintf(fp, "\n");
            else fprintf(fp, ";");
            
        }
        
    }

    // close file
	fclose(fp);
    
    // --------------------------------------------------
    
    // open file
	fp = fopen("output/output_FULL_qrand.txt", "w");
    if (!fp) {
        // print error if we cannot open file
		fprintf(stderr, "can't open file\n");
		return;
	}
    
    // loop over experiments
    for (i = 0; i < NEXP; i++) {
        
        // loop over inference times
        for (j = 0; j < g.ntinf; j++) {
            
            // write results to file
            fprintf(fp, "%u", out->qrand[i][j]);
            
            // either ";" or newline after a value
            if (j == (g.ntinf - 1)) fprintf(fp, "\n");
            else fprintf(fp, ";");
            
        }
        
    }

    // close file
	fclose(fp);
    
    // --------------------------------------------------
    
    // open file
	fp = fopen("output/output_FULL_ranks.txt", "w");
    if (!fp) {
        // print error if we cannot open file
		fprintf(stderr, "can't open file\n");
		return;
	}
    
    // loop over experiments
    for (i = 0; i < NEXP; i++) {
        
        // loop over inference times
        for (j = 0; j < g.ntinf; j++) {
            
            // write results to file
            fprintf(fp, "%u", out->ranks[i][j]);
            
            // either ";" or newline after a value
            if (j == (g.ntinf - 1)) fprintf(fp, "\n");
            else fprintf(fp, ";");
            
        }
        
    }

    // close file
	fclose(fp);
    
    // --------------------------------------------------
    
    // open file
	fp = fopen("output/output_FULL_nsources.txt", "w");
    if (!fp) {
        // print error if we cannot open file
		fprintf(stderr, "can't open file\n");
		return;
	}
    
    // loop over experiments
    for (i = 0; i < NEXP; i++) {
        
        // loop over inference times
        for (j = 0; j < g.ntinf; j++) {
            
            // write results to file
            fprintf(fp, "%u", out->nsources[i][j]);
            
            // either ";" or newline after a value
            if (j == (g.ntinf - 1)) fprintf(fp, "\n");
            else fprintf(fp, ";");
            
        }
        
    }

    // close file
	fclose(fp);
    
    // --------------------------------------------------
    
    // open file
	fp = fopen("output/output_FULL_nactive.txt", "w");
    if (!fp) {
        // print error if we cannot open file
		fprintf(stderr, "can't open file\n");
		return;
	}
    
    // loop over experiments
    for (i = 0; i < NEXP; i++) {
        
        // loop over inference times
        for (j = 0; j < g.ntinf; j++) {
            
            // write results to file
            fprintf(fp, "%u", out->nactive[i][j]);
            
            // either ";" or newline after a value
            if (j == (g.ntinf - 1)) fprintf(fp, "\n");
            else fprintf(fp, ";");
            
        }
        
    }

    // close file
	fclose(fp);
    
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// write results from QUERYING to file

void export_output_querying (OUTPUTQ *out, char s[]) {
    
    // declare unsigned integers
    unsigned int i, j, k;
    
    // declare pointer to FILE
	FILE *fp;
    
    // initialize string for filename
    char fname1[30];
    
    // cat string and integer
    snprintf(fname1, 30, "output/ranks_%s.txt", s);
    
    // open file
	fp = fopen(fname1, "w");
    if (!fp) {
        // print error if we cannot open file
		fprintf(stderr, "can't open file\n");
		return;
	}
    
    // loop over experiments
    for (i = 0; i < NEXP; i++) {
        
        // loop over inference times
        for (j = 0; j < g.ntinf; j++) {
            
            // write results to file (for different queries
            for (k = 0; k < NQUERIES; k++) fprintf(fp, "%u;", out->ranks[i][j][k]);
            
            // newline after last inference time
            if (j == (g.ntinf - 1)) fprintf(fp, "\n");
            
        }
        
    }

    // close file
	fclose(fp);
    
    // --------------------------------------------------
    
    // initialize string for filename
    char fname2[30];
    
    // cat string and integer
    snprintf(fname2, 30, "output/nsources_%s.txt", s);
    
    // open file
	fp = fopen(fname2, "w");
    if (!fp) {
        // print error if we cannot open file
		fprintf(stderr, "can't open file\n");
		return;
	}
    
    // loop over experiments
    for (i = 0; i < NEXP; i++) {
        
        // loop over inference times
        for (j = 0; j < g.ntinf; j++) {
            
            // write results to file (for different queries
            for (k = 0; k < NQUERIES; k++) fprintf(fp, "%u;", out->nsources[i][j][k]);
            
            // newline after last inference time
            if (j == (g.ntinf - 1)) fprintf(fp, "\n");
            
        }
        
    }

    // close file
	fclose(fp);
    
    // --------------------------------------------------
    
    // initialize string for filename
    char fname3[30];
    
    // cat string and integer
    snprintf(fname3, 30, "output/qmap_%s.txt", s);
    
    // open file
	fp = fopen(fname3, "w");
    if (!fp) {
        // print error if we cannot open file
		fprintf(stderr, "can't open file\n");
		return;
	}
    
    // loop over experiments
    for (i = 0; i < NEXP; i++) {
        
        // loop over inference times
        for (j = 0; j < g.ntinf; j++) {
            
            // write results to file (for different queries
            for (k = 0; k < NQUERIES; k++) fprintf(fp, "%u;", out->qmap[i][j][k]);
            
            // newline after last inference time
            if (j == (g.ntinf - 1)) fprintf(fp, "\n");
            
        }
        
    }

    // close file
	fclose(fp);
    
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -