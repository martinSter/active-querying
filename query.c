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
        printf("Node %d randomly sampled\n", max_node);
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

unsigned int maxp (INF *inf, PRIOR *pr, EVID *evid) {
    
    // declare integers
    unsigned int i, ns = 0, max_node;
    
    // create an array of doubles on stack
    double temp[g.n];
    
    // avoid warning
    UNUSED(pr);
    
    // count number of still possible sources
    for (i = 0; i < inf->n_sources; i++) if (inf->marginal_sources[i] != 0.0) ns++;
    
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
    printf("Node %d selected (probability: %.4f)\n", max_node, temp[max_node]);
    
    // return result
    return max_node;
    
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// UCTY

unsigned int ucty (INF *inf, PRIOR *pr, EVID *evid) {
    
    // declare integers
    unsigned int i, ns = 0, max_node;
    
    // declare doubles
    double ves, vei, ver;
    
    // create an array of doubles on stack
    double temp[g.n];
    
    // avoid warning
    UNUSED(pr);
    
    // count number of still possible sources
    for (i = 0; i < inf->n_sources; i++) if (inf->marginal_sources[i] != 0.0) ns++;
    
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
    printf("Node %d selected (vote entropy: %.4f)\n", max_node, temp[max_node]);
    
    // return result
    return max_node;
    
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// AKLD

unsigned int akld (INF *inf, PRIOR *pr, EVID *evid) {
    
    // declare integers
    unsigned int i, j, h, ns = 0, max_node;
    
    // declare and initialize (long) doubles
    double susc, infe, reco, temp1, temp2, temp3, probS, avg_div;
    
    // create an array of doubles on stack
    double temp[g.n];
    
    // count number of still possible sources
    for (i = 0; i < inf->n_sources; i++) if (inf->marginal_sources[i] != 0.0) ns++;
    
    // return NONE if there is only one possible source left
    if (ns == 1.0) {
        printf("Only one possible source left, returning NONE\n");
        return NONE;
    }
    
    // ============================================================
    
    // loop over all nodes in graph
    for (h = 0; h < g.n; h++) {
        
        // initialize avg_div to 0.0
        avg_div = 0.0;
        
        // compute logs of weighted predictions
        temp1 = (n[h].weighted_pred.s > 0.0) ? log(n[h].weighted_pred.s) : 0.0;
        temp2 = (n[h].weighted_pred.i > 0.0) ? log(n[h].weighted_pred.i) : 0.0;
        temp3 = (n[h].weighted_pred.r > 0.0) ? log(n[h].weighted_pred.r) : 0.0;
        
        // loop over sources
        for (i = 0; i < inf->n_sources; i++) {
        
            // loop over possible T's
            for (j = 0; j < pr->n_durs; j++) {
            
                // initialize to 0.0
                susc = 0.0; infe = 0.0; reco = 0.0;
                
                // compute probability of being susceptible
                probS = 1 - p[i].inf[h][j] - p[i].rec[h][j];
                
                // compute divergences for three possible outcomes (only if indiv. pred. larger than 0)
                if (probS > 0.0) susc = probS * (log(probS) - temp1);
                if (p[i].inf[h][j] > 0.0) infe = p[i].inf[h][j] * (log(p[i].inf[h][j]) - temp2);
                if (p[i].rec[h][j] > 0.0) reco = p[i].rec[h][j] * (log(p[i].rec[h][j]) - temp3);
                
                // multiply with posterior
                avg_div += ((susc + infe + reco) * inf->log_liks[i][j]);
                //avg_div += (sus + inf + rec);
            
            }
            
        }
        
        // store avg. KL-divergence in temp
        temp[h] = avg_div;
        
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
    printf("Node %d selected (avg. KL-divergence: %.4f)\n", max_node, temp[max_node]);
    
    // return result
    return max_node;
    
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// ONE-HOP

unsigned int onehop (INF *inf, PRIOR *pr, EVID *evid) {
    
    // declare integers
    unsigned int i, j, k, max_node, ns = 0;
    
    // create arrays
    double temp[g.n];
    
    // initialize all values of temp to 0.0
    memset(temp, 0.0, g.n * sizeof(double));
    
    // count number of still possible sources
    for (i = 0; i < inf->n_sources; i++) if (inf->marginal_sources[i] != 0.0) ns++;
    
    // return NONE if there is only one possible source left
    if (ns == 1.0) {
        printf("Only one possible source left, returning NONE\n");
        return NONE;
    }
    
    // ============================================================
    
    // 1) loop over all nodes in evidence to find nodes that can be contacted by evidence nodes
    // ==> Nodes that are contacted by evidence nodes before inference time.
    for (i = 0; i < evid->n_evid; i++) {
        
        // if node in evidence is infectious or recovered ...
        if (n[evid->nodes[i]].true_state != 0) {
            
            // loop over node's neighbors
            for (j = 0; j < n[evid->nodes[i]].deg; j++) {
                
                // loop over contact times with neighbor j
                for (k = 0; k < n[evid->nodes[i]].nc[j]; k++) {
                    
                    // if contact time is larger or equal to first possible starting time and smaller or equal to inference time...
                    if (n[evid->nodes[i]].t[j][k] >= (g.t_now - pr->n_durs) && n[evid->nodes[i]].t[j][k] <= g.t_now) {
                        
                        // add 1 to temp at position of neighbors
                        temp[n[evid->nodes[i]].nb[j]]++;
                        
                        // we can break the loop over contact times and go to next neighbor
                        break;
                        
                    }
                    
                }
                
            }
            
        }
        
    }
    
    // 2) loop over all possible sources in order to find nodes that have a contact to one of the nodes in evidence
    // note: it is enough to loop over possible sources since only these nodes have out-edges
    // ==> Nodes that contact evidence nodes before inference time.
    // This step is technically not necessary for undirected networks.
    for (i = 0; i < g.nps; i++) {
        
        // loop over node's neighbors
        for (j = 0; j < n[g.ps[i]].deg; j++) {
            
            // check if neighbor is in evidence
            if (n[n[g.ps[i]].nb[j]].e == 1) {
                
                // loop over contact times with neighbor j
                for (k = 0; k < n[g.ps[i]].nc[j]; k++) {
                    
                    // if contact time is larger or equal to first possible starting time and smaller or equal to inference time...
                    if (n[g.ps[i]].t[j][k] >= (g.t_now - pr->n_durs) && n[g.ps[i]].t[j][k] <= g.t_now) {
                        
                        // add 1 to temp at position of neighbors
                        temp[g.ps[i]]++;
                        
                        // we can break the loop over contact times and go to next neighbor
                        break;
                        
                    }
                    
                }
                
            }
            
        }
        
    }
    
    // set nodes in evidence to zero
    for (i = 0; i < evid->n_evid; i++) temp[evid->nodes[i]] = 0.0;
    
    // find node with max. number of contacts to activated nodes
    max_node = find_max(temp);
    
    // if function returns NONE, we sample a node randomly
    if (max_node == NONE) {
        // randomly sample node from all nodes that are not in evidence yet
        max_node = random_sampling(evid);
        // return result
        return max_node;
    }
    
    // print result to console
    printf("Node %d selected (number of contacts: %.0f)\n", max_node, temp[max_node]);
    
    // return result
    return max_node;
    
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// RANDOM

unsigned int randq (INF *inf, PRIOR *pr, EVID *evid) {
    
    // declare integers
    unsigned int i, j, max_node, ns = 0, sum_pn = 0;
    
    // create arrays
    unsigned int temp[g.n];
    unsigned int possible_nodes[g.n];
    
    // initialize all values of temp to 0
    memset(temp, 0, g.n * sizeof(unsigned int));
    
    // count number of still possible sources
    for (i = 0; i < inf->n_sources; i++) if (inf->marginal_sources[i] != 0.0) ns++;
    
    // return NONE if there is only one possible source left
    if (ns == 1.0) {
        printf("Only one possible source left, returning NONE\n");
        return NONE;
    }
    
    // ============================================================
    
    // loop over all possible sources
    for (i = 0; i < inf->n_sources; i++) {
        
        // jump to next source if current source is impossible
        if (inf->marginal_sources[i] == 0.0) continue;
        
        // loop over all nodes and increment temp by 1 if node has positive probability of being reached by possible source
        for (j = 0; j < g.n; j++) if((p[i].inf[j][pr->n_durs - 1] + p[i].rec[j][pr->n_durs - 1]) > 0.0) temp[j]++;
        
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
        printf("Node %d selected\n", max_node);
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