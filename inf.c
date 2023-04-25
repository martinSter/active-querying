// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// inference functions

// include header file
#include "run.h"

// declare external variables
extern GLOBALS g;
extern NODE *n;
extern MCS *p;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// function to perform log-sum-exp trick

void log_sum_exp (INF *inf, PRIOR *pr, EVID *evid) {
    
    // declare integers
    unsigned int i, j, h, status;
    
    // declare doubles
    double prob, max, lognorm = 0;
    
    // -----------------------------------------------------------------------------------------
    // FIRST LOOP OVER LOG_LIKS
    // compute log prior + log likelihood
    
    // loop over sources
    for (i = 0; i < inf->n_sources; i++) {
        
        // loop over possible T's
        for (j = 0; j < pr->n_durs; j++) {
            
            // initialize log_liks as log of prior probability (or INVALID if prior is 0.0)
            inf->log_liks[i][j] = (pr->priors[j] > 0.0) ? log(pr->priors[j]) : INVALID;
            
            // loop over nodes in evidence
            for (h = 0; h < evid->n_evid; h++) {

                // if log-likelihood takes on the value INVALID, break the loop over nodes
                // in that case, the current (source, T) pair is no longer possible
                if (inf->log_liks[i][j] == INVALID) break;
                
                // get status of node h
                status = n[evid->nodes[h]].true_state;
                
                // get node probability depending on status of node
                if (status == 0) prob = 1 - p[i].inf[evid->nodes[h]][j] - p[i].rec[evid->nodes[h]][j];
                else if (status == 1) prob = p[i].inf[evid->nodes[h]][j];
                else prob = p[i].rec[evid->nodes[h]][j];
                
                // if prob. is larger than 0.0, we add its log, otherwise we set it to INVALID
                if (prob > 0.0) inf->log_liks[i][j] += log(prob);
                else inf->log_liks[i][j] = INVALID;
                
            }
            
        }
        
    }
    
    // -----------------------------------------------------------------------------------------
    // SECOND LOOP OVER LOG_LIKS
    // find the maximum in log_liks and count possible (source, T) pairs
    
    // initialize max to -DBL_MAX and set h back to 0
    // we use h to count number of possible (source, T) pairs
    max = -DBL_MAX; h = 0;
    
    // loop over sources
    for (i = 0; i < inf->n_sources; i++) {
        
        // loop over possible T's
        for (j = 0; j < pr->n_durs; j++) {
            
            // jump to next iteration if INVALID
            if (inf->log_liks[i][j] == INVALID) continue;
            
            // find max. element in log_liks (true max. cannot be larger than 0.0)
            if (inf->log_liks[i][j] > max) max = inf->log_liks[i][j];
        
            // increase h by 1 to count an additional possible (source, T) pair
            h++;
            
        }
        
    }
    
    // -----------------------------------------------------------------------------------------
    // THIRD LOOP OVER LOG_LIKS
    // compute the log-normalization term
    
    // loop over sources
    for (i = 0; i < inf->n_sources; i++) {
        
        // loop over possible T's
        for (j = 0; j < pr->n_durs; j++) {
            
            // jump to next iteration if INVALID
            if (inf->log_liks[i][j] == INVALID) continue;
            
            // sum up e^(b_st - B) over all sources and T's
            lognorm += exp(inf->log_liks[i][j] - max);
            
        }
        
    }
    
    // part of log-sum-exp trick
    lognorm = log(lognorm) + max;
    
    // -----------------------------------------------------------------------------------------
    // FOURTH LOOP OVER LOG_LIKS
    // compute posterior probabilities
    
    // loop over sources
    for (i = 0; i < inf->n_sources; i++) {
        
        // loop over possible T's
        for (j = 0; j < pr->n_durs; j++) {
            
            // if log_liks are INVALID...
            if (inf->log_liks[i][j] == INVALID) {
                
                // then set posterior probability to 0.0
                inf->log_liks[i][j] = 0.0;

            } else {
                
                // else, compute the rest of log-sum-exp trick to get posterior probability
                inf->log_liks[i][j] = exp(inf->log_liks[i][j] - lognorm);
                
            }
            
        }
        
    }
    
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// function to marginalize posterior distribution

void marginalize (INF *inf, PRIOR *pr) {
    
    // declare integers
    unsigned int i, j;
    
    // allocate memory to arrays for marginal posterior probabilities
    inf->marginal_sources = calloc(inf->n_sources, sizeof(double));
    inf->marginal_T = calloc(pr->n_durs, sizeof(double));
    
    // loop over sources
    for (i = 0; i < inf->n_sources; i++) {
        
        // loop over possible T's
        for (j = 0; j < pr->n_durs; j++) {
            
            // add up posterior probabilities over T's
            inf->marginal_sources[i] += inf->log_liks[i][j];
            
            // add up posterior probabilities over sources
            inf->marginal_T[j] += inf->log_liks[i][j];
            
        }
        
    }
    
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// function to compute weighted node predictions

void weighted_node_predictions (INF *inf, PRIOR *pr) {
    
    // declare integers
    unsigned int i, j, h;
    
    // initialize all weighted predictions to 0.0
    for (h = 0; h < g.n; h++) {
        
        n[h].weighted_pred.s = 0.0;
        n[h].weighted_pred.i = 0.0;
        n[h].weighted_pred.r = 0.0;
    
    }
    
    // loop over all nodes in graph
    for (h = 0; h < g.n; h++) {
        
        // loop over sources
        for (i = 0; i < inf->n_sources; i++) {
        
            // loop over possible T's
            for (j = 0; j < pr->n_durs; j++) {
                
                // jump to next (source, T) pair if posterior is 0.0
                if (inf->log_liks[i][j] == 0.0) continue;
                
                // add up individual predictions using posterior as weight
                n[h].weighted_pred.s += ((1 - p[i].inf[h][j] - p[i].rec[h][j]) * inf->log_liks[i][j]);
                n[h].weighted_pred.i += (p[i].inf[h][j] * inf->log_liks[i][j]);
                n[h].weighted_pred.r += (p[i].rec[h][j] * inf->log_liks[i][j]);
                
            }
            
        }
        
    }
    
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// function to perform full inference (on all node states)

void full_inference (unsigned int nexp, INF *inf, PRIOR *pr, OUTF *outf) {
    
    // declare integer
    unsigned int i;
    
    // declare evidence struct
    EVID evid;
    
    // set number of nodes in evidence to g.n (remember that we observe all nodes here)
    evid.n_evid = g.n;
    
    // allocate memory to evid.nodes
    evid.nodes = calloc(evid.n_evid, sizeof(unsigned int));
    
    // add nodes to evid.nodes (here simply all nodes in g.n)
    for (i = 0; i < evid.n_evid; i++) evid.nodes[i] = i;
    
    // allocate memory to log_liks array of arrays (outer array)
    inf->log_liks = malloc(inf->n_sources * sizeof(double *));
    
    // allocate memory to inner arrays
    for (i = 0; i < inf->n_sources; i++) inf->log_liks[i] = calloc(pr->n_durs, sizeof(double));

    // run log-sum-exp trick to compute posterior distribution over (source, T) pairs
    log_sum_exp(inf, pr, &evid);
    
    // print (final) posterior to file
    // export_posterior(nexp, g.n, inf, pr);
    
    // compute marginal posterior distributions
    marginalize(inf, pr);
    
    // run eval_full() to compute all evaluation measures and store them in outf
    eval_full(nexp, inf, outf);
    
    // store number of affected nodes and true source in outf
    outf->nactives[nexp] = g.ns_gt;
    outf->true_sources[nexp] = g.true_source;
    
    // free memory of inner arrays
    for (i = 0; i < inf->n_sources; i++) free(inf->log_liks[i]);
    
    // free memory of outer array, marginal posterior arrays, and nodes array in evid
    free(inf->log_liks); free(inf->marginal_sources); free(inf->marginal_T); free(evid.nodes);
    
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// function to run inference and to call querying for one strategy

void infer_query (unsigned int nexp, INF *inf, PRIOR *pr, unsigned int (*fn)(INF *inf, PRIOR *pr, EVID *evid), OUTQ *o, OUTT *ot, unsigned int print) {
    
    // declare integers
    unsigned int i, next_node;
    
    // declare start and end for storing CPU run times
    clock_t t0, t1;
    
    // record start time t0
    t0 = clock();
    
    // allocate memory to log_liks array of arrays (outer array)
    inf->log_liks = malloc(inf->n_sources * sizeof(double *));
    
    // allocate memory to inner arrays
    for (i = 0; i < inf->n_sources; i++) inf->log_liks[i] = calloc(pr->n_durs, sizeof(double));
    
    // declare evidence struct
    EVID evid;
    
    // initialize evidence with first observed node
    initialize_evidence(&evid);
    
    // loop over NQUERIES
    for (i = 0; i < NQUERIES; i++) {
        
        // run log-sum-exp trick to compute posterior distribution over (source, T) pairs (based on current evidence)
        log_sum_exp(inf, pr, &evid);
        
        // print posterior to file if PRINT == 1 (set in run.h)
        if (print == 1) export_posterior(nexp, i, inf, pr);
        
        // compute marginal posterior distributions
        marginalize(inf, pr);
        
        // compute weighted node predictions
        weighted_node_predictions(inf, pr);
    
        // run eval_query() to compute all evaluation measures and store them in o
        eval_query(nexp, i, inf, o);
        
        // query next node
        next_node = fn(inf, pr, &evid);
        
        // add queried node to evidence if queried node is not NONE
        if (next_node != NONE) update_evidence(next_node, &evid);
        
    }
    
    // free memory of inner arrays
    for (i = 0; i < inf->n_sources; i++) free(inf->log_liks[i]);
    
    // free memory of outer array and marginal posterior arrays
    free(inf->log_liks); free(inf->marginal_sources); free(inf->marginal_T); free(evid.nodes);
    
    // record end time t1
    t1 = clock();
    
    // store CPU run time
    ot->run_time[nexp] = ((double)(t1 - t0)) / CLOCKS_PER_SEC;
    
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -