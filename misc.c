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
// giving exponential random numbers with a mean reciprocal of the recovery rate

unsigned int exptime () {
	uint32_t r = pcg_32();

	if (r == 4294967295u) return 0;

    // this basically implements the inverse probability transform for exponential distr.
	return -1 * log((r + 1) / 4294967296.0) / g.mu;
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// gets the index of you in me's adjacency list

unsigned int get_index (unsigned int me, unsigned int you) {
	
    // define int i
    unsigned int i;

    // loop over neighbors of 'me' and find index i where neighbor is 'you'
    // if we find neighbor, we return its index
	for (i = 0; i < n[me].deg; i++) if (n[me].nb[i] == you) return i;

	// what follows is procedure if neighbor is not yet in 'nb'
    
    // note: alloc is initialized with all zeros
    // note: we allocate more space than we need (double deg) so that we
    // do not have to allocate more space in every step
	if (alloc[me] <= n[me].deg) {
		// get value of alloc at index 'me'
        i = alloc[me];
        // set the value at alloc[me] to double the degree if degree is positive or 
        // to 1 otherwise (i.e. if deg == 0)
		alloc[me] = (n[me].deg > 0) ? 2 * n[me].deg : 1;
        // increase the allocated space to arrays with the value in alloc[me]
		n[me].nb = realloc(n[me].nb, alloc[me] * sizeof(unsigned int));
		n[me].nc = realloc(n[me].nc, alloc[me] * sizeof(unsigned int));
        // since t is a pointer to a pointer, we use the size of a pointer to int
        // effectively, t is an array of arrays containing the contact times with each neighbor
		n[me].t = realloc(n[me].t, alloc[me] * sizeof(unsigned int *));
		// for loop that does not start at 0 but at previous value of alloc[me], see above
        for ( ; i < alloc[me]; i++) {
            // initialize new positions in nc and nb to 0
			n[me].nc[i] = n[me].nb[i] = 0;
            // initialize new arrays in t to NULL
			n[me].t[i] = NULL;
		}
	}

    // add 'you' to 'nb' and increase 'deg' by one
	n[me].nb[n[me].deg++] = you;

    // return index of new neighbor
	return n[me].deg - 1;
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// reads the network, assumes a contact list with vertex label 0, N-1,
// ordered in time. If your network has nodes with zero contacts, make sure
// that none of them is the node with largest index.
// This assumes directed networks where the first node is always the orignating
// node and the second node is the receiving node.

void read_data (FILE *fp, unsigned int dir) {
    
    // declare integers
	unsigned int i, me, you;
    
    // assert that input 'dir' is one of the possible values
    // undirected: dir = 0, directed: dir = 1
    assert(dir == 0 || dir == 1);

    // initialize number of nodes and edges in g to 0
	g.n = 0;
    g.e = 0;

	// scan the system size
    // loop over lines in fp
    // that is why the node with largest integer should have contacts,
    // if it does not, then the system size will not be collected correctly
	while (2 == fscanf(fp, "%u %u %*u\n", &me, &you)) {
		// set g.n always to largest value
        if (g.n < me) g.n = me;
		if (g.n < you) g.n = you;
        // increment number of edges by 1
        g.e++;
	}

    // adds 1 to n since node indices start at 0
	g.n++;
    
    // prints out number of nodes and number of edges
    printf("Number of nodes: %d\nNumber of edges: %d\n", g.n, g.e);

    // initialize array of NODE structs in heap
    // with calloc, we initialize integers to 0
	n = calloc(g.n, sizeof(NODE));
    
    // initialize array of size system (number of nodes)
    // with calloc we initialize integers to 0
    // this is used to manage memory allocation of nb, nc, and t
	alloc = calloc(g.n, sizeof(unsigned int));

    // sets the file position to the beginning
	rewind(fp);

	// scan the degrees (here we implicitly assume undirected network)
    // loop over lines in fp
	while (2 == fscanf(fp, "%u %u %*u\n", &me, &you)) {
		// find index of 'you' in adjacency list of 'me'
        // accounting for deg and nb is taken care of in 'get_index()'
        i = get_index(me, you);
		// increase number of contacts with node 'you' by 1
        n[me].nc[i]++;
        // if network is undirected, we need to do the same but the other way around
        if (dir == 0) {
            // find index of 'me' in adjacency list of 'you'
            // accounting for deg and nb is taken care of in 'get_index()'
            i = get_index(you, me);
            // increase number of contacts with node 'me' by 1
            n[you].nc[i]++;
        }        
	}

    // sets the file position to the beginning
	rewind(fp);

    // prepare arrays of contact times
    // loop over int me (over all nodes in network)
	for (me = 0; me < g.n; me++) {
        // loop over all contacts 'i' of node 'me'
		for (i = 0; i < n[me].deg; i++) {
            // allocate memory for array of contact times
			n[me].t[i] = malloc(n[me].nc[i] * sizeof(unsigned int));
            // set number of contacts back to 0
            // this is needed because we use nc to increment in next code block
			n[me].nc[i] = 0;
		}
	}

	// scan the times
    // loop over lines in fp (this time considering contact times)
    // since input edgelist is sorted by date, the last g.dur will be length of period
	while (3 == fscanf(fp, "%u %u %u\n", &me, &you, &g.dur)) {
		// find index of 'you' in adjacency list of 'me'
        i = get_index(me, you);
        // set contact times
		n[me].t[i][n[me].nc[i]++] = g.dur;
        // if network is undirected, we need to do the same but the other way around
        if (dir == 0) {
            // find index of 'me' in adjacency list of 'you'
            i = get_index(you, me);
            // set contact times
            n[you].t[i][n[you].nc[i]++] = g.dur;
        }
	}
    
    // print g.dur to console
    printf("Duration: %d\n", g.dur);

	// reallocate adjacency lists
    // this is needed because we might have allocated too much memory in the 'get_index()' part
    // here, we therefore allocate the correct amount of memory
    // loop over all nodes in network
	for (i = 0; i < g.n; i++) {
		n[i].nb = realloc(n[i].nb, n[i].deg * sizeof(unsigned int));
		n[i].nc = realloc(n[i].nc, n[i].deg * sizeof(unsigned int));
		n[i].t = realloc(n[i].t, n[i].deg * sizeof(unsigned int *));
        // for every node i, sort t in decreasing order of its last element
		quick(i);
	}

    // free memory
	free(alloc);
    
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// find all nodes that could be a source (out-edges)

void find_possible_sources () {
    
    // declare unsigned integer
    unsigned int i;
    
    // initialize number of possible sources (nps) to 0
    g.nps = 0;
    
    // allocate memory to g.ps
    g.ps = calloc(g.n, sizeof(unsigned int));
    
    // find all nodes with out-edges
    for (i = 0; i < g.n; i++) if (n[i].deg > 0) g.ps[g.nps++] = i;
    
    // reallocate memory
    g.ps = realloc(g.ps, g.nps * sizeof(unsigned int));
    
    // print number of nodes with out-edges to console
    printf("Number of nodes with out-edges: %d\n", g.nps);
    
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// initialize the prior on possible values of T
// type = 0 (geometric), type = 1 (delta)

void init_prior (PRIOR *pr, unsigned int type) {
    
    // declare unsigned integer
    unsigned int i;
    
    // get probability corresponding to avg. T
    double prob = 1.0 / (1.0 + g.T);
    
    // assert that input 'type' is one of the possible values
    assert(type == 0 || type == 1);
    
    // determine number of possible values for T
    pr->n_durs = NUMT * g.T + 1;
    
    // allocate memory
    pr->priors = calloc(pr->n_durs, sizeof(double));
    
    // depending on type, set values of prior
    if (type == 0) {
        
        // truncated (geometric) prior probabilities
        // CDF is 1.0 - (1.0 - p)^(cutoff + 1), hence we can plug-in pr->n_durs (no need to subtract 1)
        for (i = 0; i < pr->n_durs; i++) pr->priors[i] = (prob * pow(1.0 - prob, i)) / (1.0 - pow(1.0 - prob, pr->n_durs));
        
    } else if (type == 1) {
        
        // Delta prior (only 1.0 for true T)
        for (i = 0; i < pr->n_durs; i++) if (i == g.T) pr->priors[i] = 1.0;
        
    }
    
}


    
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// set up output structs for querying

void set_up (OUTQ *out) {
    
    // declare integer
    unsigned int i;
    
    // allocate memory to arrays in out
    out->ranks = malloc(NEXP * sizeof(unsigned int *));
    out->nsources = malloc(NEXP * sizeof(unsigned int *));
    out->css = malloc(NEXP * sizeof(unsigned int *));
    
    // loop over number of experiments
    for (i = 0; i < NEXP; i++) {
        
        // allocate memory to arrays on inner level
        out->ranks[i] = malloc(NQUERIES * sizeof(unsigned int));
        out->nsources[i] = malloc(NQUERIES * sizeof(unsigned int));
        out->css[i] = malloc(NQUERIES * sizeof(unsigned int));
        
    }
    
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// clean up output structs for querying

void clean_up (OUTQ *out) {
    
    // declare integer
    unsigned int i;
    
    // loop over number of experiments
    for (i = 0; i < NEXP; i++) {
        
        // free arrays on inner level
        free(out->ranks[i]);
        free(out->nsources[i]);
        free(out->css[i]);
        
    }
    
    // free arrays on outer level
    free(out->ranks); free(out->nsources); free(out->css);
    
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// drawing the first observed infected node

void draw_first_node () {
    
    // declare integers
    unsigned int i, active = 0;
    
    // create array for active nodes
    unsigned int active_nodes[g.n];
    
    // loop over all nodes in graph and add them to array if they have been activated
    for (i = 0; i < g.n; i++) if (n[i].true_state == 1 || n[i].true_state == 2) active_nodes[active++] = i;

    // random draw of integer in active nodes
    i = pcg_32_bounded(active);
    
    // store first observed
    g.fin = active_nodes[i];
    
    // print to console
    printf("First observed node is %u\n", g.fin);
    
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
// compute measure to evaluate performance (FULL INFERENCE)

void eval_full (unsigned int nexp, INF *inf, OUTF *outf) {
    
    // declare integers
    unsigned int i, j, temp1;
    
    // declare doubles
    double temp2;
    
    // initialize arrays on the stack
    unsigned int source_indices[inf->n_sources]; // this is needed for the sorting
    double post[inf->n_sources]; // this is needed for the sorting
    double ranks[inf->n_sources]; // ranks can be decimal numbers
    
    // --------------------------------------------------
    // SORT MARGINAL SOURCE POSTERIOR IN DESCENDING ORDER
    
    // loop over sources to initialize sorting
    for (i = 0; i < inf->n_sources; i++) {
        
        // fill source_indices with i
        source_indices[i] = i;
        
        // fill post with posterior probabilities
        post[i] = inf->marginal_sources[i];
        
    }
 
    // loop over array elements
    for (i = 0; i < inf->n_sources; i++) {
        
        // loop over all array elements after i-th element
        for (j = i + 1; j < inf->n_sources; j++) {
            
            // swap elements if i-th element is smaller than j-th element
            if (post[i] < post[j]) {
                
                // store source and posterior at i in temp
                temp1 = source_indices[i];
                temp2 = post[i];
                
                // set i to j
                source_indices[i] = source_indices[j];
                post[i] = post[j];
                
                // set j from temp
                source_indices[j] = temp1;
                post[j] = temp2;
                
            }
        }
    }
    
    // -------------------------------------------
    // LOOP OVER SORTED POSTERIOR TO COMPUTE RANKS
    
    // initialize rank to 1
    unsigned int rank = 1, nn;
    
    // set i back to 0
    i = 0;
    
    // loop over array
    while (i < inf->n_sources) {
        // set j to current value of i
        j = i;
        // count how many elements with same log-liks
        while (j < inf->n_sources - 1 && post[j] == post[j + 1]) j++;
        // set nn
        nn = j - i + 1;
        // set ranks
        for (j = 0; j < nn; j++) ranks[j + i] = rank + (nn - 1) * 0.5;
        // increment rank and i
        rank += nn;
        i += nn;
    }

    // -------------------------------------------
    // LOOP OVER SORTED POSTERIOR
    
    // initialize credible set size (css) and number of active sources (non-zero posterior) to 0
    unsigned int css = 0, nact = 0;
    
    // initialize cumulated posterior to 0.0
    double cumpost = 0.0;
    
    // loop over sorted posterior
    for (i = 0; i < inf->n_sources; i++) {
        
        // increment nact by 1 as long as posterior is non-zero
        if (post[i] > 0.0) nact++;
        
        // check if cumpost is smaller than 0.95
        if (cumpost < 0.95) {
            
            // add current posterior prob. to cumulated posterior
            cumpost += post[i];
            
            // increment credible set size by 1
            css++;
            
        }
        
        // check if source at current index is true source
        if (inf->sources[source_indices[i]] == g.true_source) {
            
            // add rank of true source to outf
            outf->ranks[nexp] = ranks[i];
            
            // add posterior of true source to outf
            outf->posteriors[nexp] = post[i];
            
        }

    }
    
    // add number of active sources and credible set size to outf
    outf->nsources[nexp] = nact; outf->css[nexp] = css;
    
    // print info to console
    printf("Rank of true source: %d (Credible set size: %d)\n", (int) outf->ranks[nexp], css);
    
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// compute measure to evaluate performance (QUERYING)

void eval_query (unsigned int nexp, unsigned int nquery, INF *inf, OUTQ *o) {
    
    // declare integers
    unsigned int i, j, temp1;
    
    // declare doubles
    double temp2;
    
    // initialize arrays on the stack
    unsigned int source_indices[inf->n_sources]; // this is needed for the sorting
    double post[inf->n_sources]; // this is needed for the sorting
    double ranks[inf->n_sources]; // ranks can be decimal numbers
    
    // --------------------------------------------------
    // SORT MARGINAL SOURCE POSTERIOR IN DESCENDING ORDER
    
    // loop over sources to initialize sorting
    for (i = 0; i < inf->n_sources; i++) {
        
        // fill source_indices with i
        source_indices[i] = i;
        
        // fill post with posterior probabilities
        post[i] = inf->marginal_sources[i];
        
    }
 
    // loop over array elements
    for (i = 0; i < inf->n_sources; i++) {
        
        // loop over all array elements after i-th element
        for (j = i + 1; j < inf->n_sources; j++) {
            
            // swap elements if i-th element is smaller than j-th element
            if (post[i] < post[j]) {
                
                // store source and posterior at i in temp
                temp1 = source_indices[i];
                temp2 = post[i];
                
                // set i to j
                source_indices[i] = source_indices[j];
                post[i] = post[j];
                
                // set j from temp
                source_indices[j] = temp1;
                post[j] = temp2;
                
            }
        }
    }
    
    // -------------------------------------------
    // LOOP OVER SORTED POSTERIOR TO COMPUTE RANKS
    
    // initialize rank to 1
    unsigned int rank = 1, nn;
    
    // set i back to 0
    i = 0;
    
    // loop over array
    while (i < inf->n_sources) {
        // set j to current value of i
        j = i;
        // count how many elements with same log-liks
        while (j < inf->n_sources - 1 && post[j] == post[j + 1]) j++;
        // set nn
        nn = j - i + 1;
        // set ranks
        for (j = 0; j < nn; j++) ranks[j + i] = rank + (nn - 1) * 0.5;
        // increment rank and i
        rank += nn;
        i += nn;
    }

    // -------------------------------------------
    // LOOP OVER SORTED POSTERIOR
    
    // initialize credible set size (css) and number of active sources (non-zero posterior) to 0
    unsigned int css = 0, nact = 0;
    
    // initialize cumulated posterior to 0.0
    double cumpost = 0.0;
    
    // loop over sorted posterior
    for (i = 0; i < inf->n_sources; i++) {
        
        // increment nact by 1 as long as posterior is non-zero
        if (post[i] > 0.0) nact++;
        
        // check if cumpost is smaller than 0.95
        if (cumpost < 0.95) {
            
            // add current posterior prob. to cumulated posterior
            cumpost += post[i];
            
            // increment credible set size by 1
            css++;
            
        }
        
        // check if source at current index is true source; if yes, add rank to o
        if (inf->sources[source_indices[i]] == g.true_source) o->ranks[nexp][nquery] = ranks[i];

    }
    
    // add number of active sources and credible set size to o
    o->nsources[nexp][nquery] = nact; o->css[nexp][nquery] = css;
    
    // print info to console
    printf("Rank of true source: %d (Credible set size: %d)\n", (int) o->ranks[nexp][nquery], css);
    
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// write posterior to file

void export_posterior (unsigned int nexp, unsigned int nq, INF *inf, PRIOR *pr) {
    
    // declare unsigned integers
    unsigned int i, j;
    
    // declare pointer to FILE
	FILE *fp;
    
    // initialize string for filename
    char fname[30];
    
    // cat string and integer
    snprintf(fname, 30, "POSTERIOR_EXP_%u_Q_%u.txt", nexp, nq);
    
    // open file
	fp = fopen(fname, "w");
    if (!fp) {
        // print error if we cannot open file
		fprintf(stderr, "can't open file\n");
		return;
	}
    
    // add empty because first column will contain sources
    fprintf(fp, ";");
    
    // print T's in top row
    for (j = 0; j < pr->n_durs; j++) fprintf(fp, "%u;", j);
    
    // add line break after top row
    fprintf(fp, "\n");
    
    // loop over sources
    for (i = 0; i < inf->n_sources; i++) {
        
        // add source as first value in row
        fprintf(fp, "%u;", inf->sources[i]);
        
        // loop over possible T's
        for (j = 0; j < pr->n_durs; j++) {
            
            // add posterior probability
            fprintf(fp, "%f;", inf->log_liks[i][j]);
            
        }
        
        // add line break at the end
        fprintf(fp, "\n");
        
    }

    // close file
	fclose(fp);
    
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// write results to file

void export_output (OUTF *outf) {
    
    // declare unsigned integers
    unsigned int i;
    
    // declare pointer to FILE
	FILE *fp;
    
    // initialize string for filename
    char fname1[15];
    
    // cat string and integer
    snprintf(fname1, 15, "FULL_%u.txt", g.T);
    
    // open file
	fp = fopen(fname1, "w");
    if (!fp) {
        // print error if we cannot open file
		fprintf(stderr, "can't open file\n");
		return;
	}
    
    // loop over experiments and write results to file
    for (i = 0; i < NEXP; i++) fprintf(fp, "%u;%u;%u;%u;%f;%u;%f\n", outf->nactives[i], outf->nsources[i], outf->true_sources[i], outf->ranks[i], outf->posteriors[i], outf->css[i], outf->sim_time[i]);

    // close file
	fclose(fp);
    
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// write results from QUERYING to file

void export_output_querying (OUTQ *o, char s[]) {
    
    // declare unsigned integers
    unsigned int i, j;
    
    // declare pointer to FILE
	FILE *fp;
    
    // initialize string for filename
    char fname1[20];
    
    // cat string and integer
    snprintf(fname1, 20, "%s_ranks.txt", s);
    
    // open file
	fp = fopen(fname1, "w");
    if (!fp) {
        // print error if we cannot open file
		fprintf(stderr, "can't open file\n");
		return;
	}
    
    // loop over experiments
    for (i = 0; i < NEXP; i++) {
        
        // loop over queries and write ranks to file
        for (j = 0; j < NQUERIES; j++) fprintf(fp, "%u;", o->ranks[i][j]);
        
        // print new line before printing results of next experiment
        fprintf(fp, "\n");
        
    }

    // close file
	fclose(fp);
    
    // --------------------------------------------------
    
    // initialize string for filename
    char fname2[20];
    
    // cat string and integer
    snprintf(fname2, 20, "%s_nsources.txt", s);
    
    // open file
	fp = fopen(fname2, "w");
    if (!fp) {
        // print error if we cannot open file
		fprintf(stderr, "can't open file\n");
		return;
	}
    
    // loop over experiments
    for (i = 0; i < NEXP; i++) {
        
        // loop over queries and write ranks to file
        for (j = 0; j < NQUERIES; j++) fprintf(fp, "%u;", o->nsources[i][j]);
        
        // print new line before printing results of next experiment
        fprintf(fp, "\n");
        
    }

    // close file
	fclose(fp);
    
    // --------------------------------------------------
    
    // initialize string for filename
    char fname3[20];
    
    // cat string and integer
    snprintf(fname3, 20, "%s_css.txt", s);
    
    // open file
	fp = fopen(fname3, "w");
    if (!fp) {
        // print error if we cannot open file
		fprintf(stderr, "can't open file\n");
		return;
	}
    
    // loop over experiments
    for (i = 0; i < NEXP; i++) {
        
        // loop over queries and write ranks to file
        for (j = 0; j < NQUERIES; j++) fprintf(fp, "%u;", o->css[i][j]);
        
        // print new line before printing results of next experiment
        fprintf(fp, "\n");
        
    }

    // close file
	fclose(fp);
    
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// write run times of querying strategies to file

void export_run_times(OUTT *ot_maxp, OUTT *ot_ucty, OUTT *ot_akld, OUTT *ot_onehop, OUTT *ot_randq) {
    
    // declare unsigned integers
    unsigned int i;
    
    // declare pointer to FILE
	FILE *fp;
    
    // initialize string for filename
    char fname[18];
    
    // cat string and integer
    snprintf(fname, 18, "RUN_TIMES_%u.txt", g.T);
    
    // open file
	fp = fopen(fname, "w");
    if (!fp) {
        // print error if we cannot open file
		fprintf(stderr, "can't open file\n");
		return;
	}
    
    // loop over experiments and write results to file
    for (i = 0; i < NEXP; i++) fprintf(fp, "%f;%f;%f;%f;%f\n", ot_maxp->run_time[i], ot_ucty->run_time[i], ot_akld->run_time[i], ot_onehop->run_time[i], ot_randq->run_time[i]);

    // close file
	fclose(fp);
    
    
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -