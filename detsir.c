// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// code for finding all possible sources

// include header file
#include "run.h"

// declare external variables
extern GLOBALS g;
extern NODE *n;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// function to find the earliest contact from now with a given neighbor

unsigned int earliest_contact (unsigned int *t, unsigned int nt, unsigned int now) {
	
    // declare some integers
    // hi = nt - 1 because index starts at 0
    unsigned int lo = 0, mid, hi = nt - 1;

    // no need to search further because t is sorted
    // return END which is max unsigned integer
	if (t[hi] <= now) return END;

	// the actual bisection search
	do {
        // >> is a right shift operator (here bits of '(lo + hi)' are shifted to the right by 1 place)
        // effectively, this computes floor((lo + hi) / 2)
		mid = (lo + hi) >> 1;
		if (t[mid] > now) hi = mid;
		else lo = mid;
	} while (hi > lo + 1);

    // the only case lo is correct
    // if now is smaller than first contact time
	if (now < t[lo]) hi = lo;

	// return the time of the contact
	return t[hi];
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// function for deterministic infection

void deterministic_infect () {
    
    // declare integers, me is first element of heap
	unsigned int i, you, t, me = g.heap[1];
    
    // declare integers, now is the infection time of me, duration is recovery time of me
	unsigned int now = n[me].time;

    // take the newly infected off the heap
	del_root();

    // go through the neighbors of the infected node
	for (i = 0; i < n[me].deg; i++) {
            
        // get neighbor i as you
		you = n[me].nb[i];
            
        // if you is S (susceptible), you can be infected
        // S is a macro defined in header file
		if (S(you)) {
                
            // find the infection time of you
            t = earliest_contact(n[me].t[i], n[me].nc[i], now);
                
            // bcoz the sorting of nbs, we can break
            // hmmm, so neighbors seem to be sorted so that we don't need to continue for loop
			if (t == END) break;
            
            // only accept infections that are at the latest at t_now
            if (t > g.t_now) continue;

			// if t is smaller than a potentially listed previous infection event, modify t
			if (t < n[you].time) {
					
                // set you's infection time
                n[you].time = t;
                    
                // if not listed before, then extend the heap
                if (n[you].heap == NONE) {
					g.heap[++g.nheap] = you;
					n[you].heap = g.nheap;
				}
                    
                // this works because there the only heap relationship that can be 
                // violated is the one between you and its parent
                // here, we upheap if necessary
				up_heap(n[you].heap);
				
            }
		}
	}

    // mark me as visited
    n[me].visited = 1;
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void deterministic_spread (unsigned int source, unsigned int start) {
    
    // declare unsigned integers
	unsigned int i;
    
    // set all nodes to not visited
    for (i = 0; i < g.n; i++) n[i].visited = 0;
    
    // infect the source at time start (0, ..., g.dur)
	n[source].time = start;
	
    // set source on position 1 of heap
    n[source].heap = 1;
    
    // add source to heap
	g.heap[g.nheap = 1] = source;

	// run the outbreak as long as heap contains elements
	while (g.nheap) deterministic_infect();

	// set time of infected and recovered nodes back to NONE
	for (i = 0; i < g.n; i++) n[i].heap = n[i].time = NONE;
    
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// find all nodes that could infect all nodes in ground truth

void find_active_sources (INF *inf, unsigned int start) {
    
    // declare some integers
    unsigned int i;
    
    // initialize number of sources to 0
    inf->n_sources = 0;
    
    // allocate memory to sources array
    inf->sources = calloc(g.nps, sizeof(unsigned int));
    
    // loop over all nodes that have out-edges
    for (i = 0; i < g.nps; i++) {
        
        // deterministic spread from source at position i
        deterministic_spread(g.ps[i], start);
        
        // if first detected node is visited, source is added as a possible source
        if (n[g.fin].visited == 1) inf->sources[inf->n_sources++] = g.ps[i];
        
    }
    
    // reallocate memory
    inf->sources = realloc(inf->sources, inf->n_sources * sizeof(unsigned int));
    
    // print to console
    printf("There are %d possible sources.\n", inf->n_sources);
    
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -