// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// active querying functions

// include header file
#include "run.h"

// declare external variables
extern GLOBALS g;
extern NODE *n;
extern MCS *p;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// random sampling

unsigned int random_sampling (EVID *evid) {
    
    // declare integers
    unsigned int i, max_node, sum_pn = 0;
    
    // create arrays
    unsigned int temp[g.n];
    unsigned int possible_nodes[g.n];
    
    // initialize all values of temp to 1
    memset(temp, 1, g.n * sizeof(unsigned int));
    
    // set nodes in evidence to zero
    for (i = 0; i < evid->n_evid; i++) temp[evid->nodes[i]] = 0;
    
    // select only nodes that are not zero
    for (i = 0; i < g.n; i++) if (temp[i] > 0) possible_nodes[sum_pn++] = i;
    
    // return node
    if (sum_pn > 0) {
        // randomly pick one
        max_node = possible_nodes[pcg_32_bounded(sum_pn)];
        // print message to console
        // printf("Node %d randomly sampled\n", max_node);
        // return node
        return max_node;
    }
    else {
        // print message to console
        printf("No more nodes available, returning NONE\n");
        return NONE;
    }
    
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// MAXP

unsigned int maxp (unsigned int tinf_idx, POSTERIOR *post, EVID *evid) {
    
    // declare integers
    unsigned int i, ns = 0, max_node;
    
    // avoid warning message
    UNUSED(tinf_idx);
    
    // create an array of doubles on stack
    double temp[g.n];
    
    // count number of still possible sources
    for (i = 0; i < g.n; i++) if (post->log_liks[i] != INVALID) ns++;
    
    // return NONE if there is only one possible source left
    if (ns == 1.0) {
        printf("Only one possible source left, returning NONE\n");
        return NONE;
    }
    
    // ============================================================
    
    // set elements in temp to probability of not being susceptible for every node
    for (i = 0; i < g.n; i++) temp[i] = 1.0 - n[i].weighted_pred.s;
    
    // set prob. of nodes in evidence to zero
    for (i = 0; i < evid->n_evid; i++) temp[evid->nodes[i]] = 0.0;
    
    // find node with max. prob.
    max_node = find_max(temp);
    
    // if function returns NONE (no node has weighted prob. > 0), we sample a node randomly
    if (max_node == NONE) {
        // randomly sample node from all nodes that are not in evidence yet
        max_node = random_sampling(evid);
        // return result
        return max_node;
    }
    
    // print result to console
    // printf("Node %d selected (probability: %.4f)\n", max_node, temp[max_node]);
    
    // return result
    return max_node;
    
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// UCTY

unsigned int ucty (unsigned int tinf_idx, POSTERIOR *post, EVID *evid) {
    
    // declare integers
    unsigned int i, ns = 0, max_node;
    
    // declare doubles
    double ves, vei, ver;
    
    // avoid warning message
    UNUSED(tinf_idx);
    
    // create an array of doubles on stack
    double temp[g.n];
    
    // count number of still possible sources
    for (i = 0; i < g.n; i++) if (post->log_liks[i] != INVALID) ns++;
    
    // return NONE if there is only one possible source left
    if (ns == 1.0) {
        printf("Only one possible source left, returning NONE\n");
        return NONE;
    }
    
    // ============================================================
    
    // loop over all nodes in graph
    for (i = 0; i < g.n; i++) {
        
        // compute individual terms
        ves = (n[i].weighted_pred.s > 0) ? n[i].weighted_pred.s * log(n[i].weighted_pred.s) : 0.0;
        vei = (n[i].weighted_pred.i > 0) ? n[i].weighted_pred.i * log(n[i].weighted_pred.i) : 0.0;
        ver = (n[i].weighted_pred.r > 0) ? n[i].weighted_pred.r * log(n[i].weighted_pred.r) : 0.0;
        
        // add up terms
        temp[i] = - (ves + vei + ver);
        
    }
    
    // set temp to 0 for all observed nodes
    for (i = 0; i < evid->n_evid; i++) temp[evid->nodes[i]] = 0.0;
    
    // find node with max. soft vote entropy
    max_node = find_max(temp);
    
    // if function returns NONE (no node has entropy > 0), we sample a node randomly
    if (max_node == NONE) {
        // randomly sample node from all nodes that are not in evidence yet
        max_node = random_sampling(evid);
        // return result
        return max_node;
    }
    
    // print result to console
    // printf("Node %d selected (vote entropy: %.4f)\n", max_node, temp[max_node]);
    
    // return result
    return max_node;
    
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// AKLD

unsigned int akld (unsigned int tinf_idx, POSTERIOR *post, EVID *evid) {
    
    // declare integers
    unsigned int i, j, ns = 0, max_node;
    
    // declare and initialize (long) doubles
    double sus, inf, rec, temp1, temp2, temp3, probS, avg_div;
    
    // create an array of doubles on stack
    double temp[g.n];
    
    // count number of still possible sources
    for (i = 0; i < g.n; i++) if (post->log_liks[i] != INVALID) ns++;
    
    // return NONE if there is only one possible source left
    if (ns == 1.0) {
        printf("Only one possible source left, returning NONE\n");
        return NONE;
    }
    
    // ============================================================
    
    // loop over all nodes in graph
    for (i = 0; i < g.n; i++) {
        
        // initialize avg_div to 0.0
        avg_div = 0.0;
        
        // compute logs of weighted predictions
        temp1 = (n[i].weighted_pred.s > 0.0) ? log(n[i].weighted_pred.s) : 0.0;
        temp2 = (n[i].weighted_pred.i > 0.0) ? log(n[i].weighted_pred.i) : 0.0;
        temp3 = (n[i].weighted_pred.r > 0.0) ? log(n[i].weighted_pred.r) : 0.0;
        
        // loop over sources (predictions)
        for (j = 0; j < g.n; j++) {
            
            // initialize to 0.0
            sus = 0.0; inf = 0.0; rec = 0.0;
            
            // compute probability of being susceptible
            probS = 1 - p[j].inf[i][tinf_idx] - p[j].rec[i][tinf_idx];
            
            // compute divergences for three possible outcomes (only if indiv. pred. larger than 0)
            if (probS > 0.0) sus = probS * (log(probS) - temp1);
            if (p[j].inf[i][tinf_idx] > 0.0) inf = p[j].inf[i][tinf_idx] * (log(p[j].inf[i][tinf_idx]) - temp2);
            if (p[j].rec[i][tinf_idx] > 0.0) rec = p[j].rec[i][tinf_idx] * (log(p[j].rec[i][tinf_idx]) - temp3);
            
            // multiply with posterior
            avg_div += ((sus + inf + rec) * post->posterior[j]);
            //avg_div += (sus + inf + rec);
            
        }
        
        // store avg. KL-divergence in temp
        temp[i] = avg_div;
        
    }
    
    // set temp to 0 for all observed nodes
    for (i = 0; i < evid->n_evid; i++) temp[evid->nodes[i]] = 0.0;
    
    // find node with max. KL-divergence
    max_node = find_max(temp);
    
    // if function returns NONE (no node has avg. KL-divergence > 0), we sample a node randomly
    if (max_node == NONE) {
        // randomly sample node from all nodes that are not in evidence yet
        max_node = random_sampling(evid);
        // return result
        return max_node;
    }
    
    // print result to console
    // printf("Node %d selected (avg. KL-divergence: %.4f)\n", max_node, temp[max_node]);
    
    // return result
    return max_node;
    
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// ONE-HOP

unsigned int onehop (unsigned int tinf_idx, POSTERIOR *post, EVID *evid) {
    
    // declare integers
    unsigned int i, j, max_node, ns = 0, sum_pn = 0;
    
    // create arrays
    unsigned int temp[g.n];
    unsigned int possible_nodes[g.n];
    
    // initialize all values of temp to 0
    memset(temp, 0, g.n * sizeof(unsigned int));
    
    // count number of still possible sources
    for (i = 0; i < g.n; i++) if (post->log_liks[i] != INVALID) ns++;
    
    // return NONE if there is only one possible source left
    if (ns == 1.0) {
        printf("Only one possible source left, returning NONE\n");
        return NONE;
    }
    
    // ============================================================
    
    // 1) loop over all nodes in evidence to find nodes that can be contacted by evidence nodes
    for (i = 0; i < evid->n_evid; i++) {
        
        // if node in evidence is infectious or recovered ...
        if (n[evid->nodes[i]].ground_truth[tinf_idx] != 0) {
            
            // loop over node's neighbors and add 1 to temp at position of neighbors
            for (j = 0; j < n[evid->nodes[i]].deg; j++) temp[n[evid->nodes[i]].nb[j]]++;
            
        }
        
    }
    
    // set nodes in evidence to zero
    for (i = 0; i < evid->n_evid; i++) temp[evid->nodes[i]] = 0;
    
    // select only nodes that are not zero
    for (i = 0; i < g.n; i++) if (temp[i] > 0) possible_nodes[sum_pn++] = i;
    
    // return node
    if (sum_pn > 0) {
        // randomly pick one
        max_node = possible_nodes[pcg_32_bounded(sum_pn)];
        // print message to console
        // printf("Node %d selected\n", max_node);
        // return node
        return max_node;
    }
    else {
        // randomly sample node from all nodes that are not in evidence yet
        max_node = random_sampling(evid);
        // return result
        return max_node;
    }
    
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// RANDOM

unsigned int randq (unsigned int tinf_idx, POSTERIOR *post, EVID *evid) {
    
    // declare integers
    unsigned int i, j, max_node, ns = 0, sum_pn = 0;
    
    // create arrays
    unsigned int temp[g.n];
    unsigned int possible_nodes[g.n];
    
    // initialize all values of temp to 0
    memset(temp, 0, g.n * sizeof(unsigned int));
    
    // count number of still possible sources
    for (i = 0; i < g.n; i++) if (post->log_liks[i] != INVALID) ns++;
    
    // return NONE if there is only one possible source left
    if (ns == 1.0) {
        printf("Only one possible source left, returning NONE\n");
        return NONE;
    }
    
    // ============================================================
    
    // loop over all possible sources
    for (i = 0; i < g.n; i++) {
        
        // jump to next source if current source is impossible
        if (post->log_liks[i] == INVALID) continue;
        
        // loop over all nodes and increment temp by 1 if node has positive probability of being reached by possible source
        for (j = 0; j < g.n; j++) if((p[i].inf[j][tinf_idx] + p[i].rec[j][tinf_idx]) > 0.0) temp[j]++;
        
    }
    
    // set nodes in evidence to zero
    for (i = 0; i < evid->n_evid; i++) temp[evid->nodes[i]] = 0;
    
    // select only nodes that are not zero
    for (i = 0; i < g.n; i++) if (temp[i] > 0) possible_nodes[sum_pn++] = i;
    
    // return node
    if (sum_pn > 0) {
        // randomly pick one
        max_node = possible_nodes[pcg_32_bounded(sum_pn)];
        // print message to console
        // printf("Node %d selected\n", max_node);
        // return node
        return max_node;
    }
    else {
        // randomly sample node from all nodes that are not in evidence yet
        max_node = random_sampling(evid);
        // return result
        return max_node;
    }
    
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -