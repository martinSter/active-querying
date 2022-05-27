// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// inference functions

// include header file
#include "run.h"

// declare external variables
extern GLOBALS g;
extern NODE *n;
extern MCS *p;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// function to perform full inference (on all node states)

void full_inference (unsigned int nexp, POSTERIOR *post, OUTPUT *o) {
    
    // declare integers
    unsigned int j, h, k, active, status;
    
    // declare double
    double prob;
    
    // loop over inference times
    for (j = 0; j < g.ntinf; j++) {
        
        // allocate memory        
        post->log_liks = calloc(g.n, sizeof(double));
        
        // initialize active to 0
        active = 0;
        
        // count number of activated nodes at inference time just
        for (h = 0; h < g.n; h++) if (n[h].ground_truth[j] > 0) active++;
        
        // print info to console
        printf("Activated nodes at time %3.2f: %d\n", g.tinf[j], active);
    
        // loop over sources
        for (k = 0; k < g.n; k++) {
            
            // loop over all nodes
            for (h = 0; h < g.n; h++) {
                
                // if log-likelihood takes on the value INVALID, break the loop over nodes
                // in that case, the current source k is no longer a possible source
                if (post->log_liks[k] == INVALID) break;
                
                // get status of node h at inference time j
                status = n[h].ground_truth[j];
                
                // get node probability depending on status of node
                if (status == 0) prob = 1 - p[k].inf[h][j] - p[k].rec[h][j];
                else if (status == 1) prob = p[k].inf[h][j];
                else prob = p[k].rec[h][j];
                
                // if prob. is larger than 0.0, we add its log, otherwise we set it to INVALID
                if (prob > 0.0) post->log_liks[k] += log(prob);
                else post->log_liks[k] = INVALID;
                
            }
            
        }
            
        // compute rank of true source
        sort_rank(post);
        
        // print result to console
        // printf("Rank of true source: %d (%d)\n", (int) g.rank, g.nsources);
        
        // store results in out_full
        o->qmap[nexp][j] = g.qmap;
        o->ranks[nexp][j] = g.rank;
        o->nsources[nexp][j] = g.nsources;
        o->nactive[nexp][j] = active;
        
        // randomly select a node from all activated nodes
        o->qrand[nexp][j] = random_draw(j);
        
        // free memory
        free(post->log_liks);
        
    }
    
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// function to compute likelihood based on evidence

void likelihood (unsigned int tinf_idx, POSTERIOR *post, EVID *evid) {
    
    // declare integers
    unsigned int i, q, status;
    
    // declare double
    double prob;
    
    // loop over sources
    for (i = 0; i < g.n; i++) {
        
        // if log-likelihood takes on the value INVALID, break the loop over nodes
        // in that case, the current source k is no longer a possible source
        if (post->log_liks[i] == INVALID) continue;
        
        // get queried node (last node in evid)
        q = evid->nodes[evid->n_evid - 1];
        
        // get status of queried node (at inference time at index tinf_idx)
        status = n[q].ground_truth[tinf_idx];
        
        // get node probability depending on status of node
        if (status == 0) prob = 1 - p[i].inf[q][tinf_idx] - p[i].rec[q][tinf_idx];
        else if (status == 1) prob = p[i].inf[q][tinf_idx];
        else prob = p[i].rec[q][tinf_idx];
        
        // if prob. is larger than 0.0, we add its log, otherwise we set it to INVALID
        if (prob > 0.0) post->log_liks[i] += log(prob);
        else post->log_liks[i] = INVALID;
        
    }
    
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// function to compute posterior

void posterior (POSTERIOR *post) {
    
    // declare integers
    unsigned int i;
    
    // declare doubles
    double thres, temp, norm = 0.0;
    
    // compute threshold
    thres = log(PREC) - log(g.nsources);
    
    // loop over sources
    for (i = 0; i < g.n; i++) {
        
        // we only need to manipulate log-likelihoods that are not INVALID
        // this is ok since posterior array is initialized with 0s (calloc)
        if (post->log_liks[i] != INVALID) {
            
            // subtract max. log-likelihood from each log value
            temp = post->log_liks[i] - g.maxlog;
        
            // exponentiate if temp is larger or equal to thres, otherwise set to 0.0
            post->posterior[i] = (temp >= thres) ? exp(temp) : 0.0;
        
        }
        
    }
    
    // sum up all likelihoods
    for (i = 0; i < g.n; i++) norm += post->posterior[i];
    
    // normalize in order to get posterior probabilities
    for (i = 0; i < g.n; i++) post->posterior[i] = post->posterior[i] / norm;
    
    // print posterior of true source to console
    // printf("Posterior of true source: %f\n", post->posterior[g.true_source]);
    
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// function to compute weighted node predictions

void weighted_node_predictions (unsigned int tinf_idx, POSTERIOR *post) {
    
    // declare integers
    unsigned int i, j;
    
    // initialize all weighted predictions to 0.0
    for (i = 0; i < g.n; i++) {
        n[i].weighted_pred.s = 0.0;
        n[i].weighted_pred.i = 0.0;
        n[i].weighted_pred.r = 0.0;
    }
    
    // loop over all nodes in graph
    for (i = 0; i < g.n; i++) {
        
        // loop over all node predictions (number of sources)
        // multiply node prediction of source with posterior of this source
        
        // loop over sources
        for (j = 0; j < g.n; j++) {
            
            // jump to next source if posterior is 0.0
            if (post->posterior[j] == 0.0) continue;
            
            // add up individual predictions and weight with posterior
            n[i].weighted_pred.s += ((1 - p[j].inf[i][tinf_idx] - p[j].rec[i][tinf_idx]) * post->posterior[j]);
            n[i].weighted_pred.i += (p[j].inf[i][tinf_idx] * post->posterior[j]);
            n[i].weighted_pred.r += (p[j].rec[i][tinf_idx] * post->posterior[j]);
            
        }
       
    }
    
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// function to run inference (at some inference time) and to call querying for one strategy

void infer_query (unsigned int nexp, unsigned int tinf_idx, POSTERIOR *post, unsigned int (*fn)(unsigned int tinf_idx, POSTERIOR *post, EVID *evid), OUTPUTQ *o) {
    
    // declare integers
    unsigned int i, next_node, trigger = 1;
    
    // declare evidence struct
    EVID evid;
    
    // allocate memory        
    post->log_liks = calloc(g.n, sizeof(double));
    post->posterior = calloc(g.n, sizeof(double));
    
    // initialize evidence with first observed node
    initialize_evidence(&evid);
    
    // loop over NQUERIES
    for (i = 0; i < NQUERIES; i++) {
        
        // compute log-likelihood (only if trigger is set to 1)
        // this is important since likelihood should not be updated if no node is queried
        if (trigger == 1) likelihood(tinf_idx, post, &evid);
        
        // compute rank of true source
        sort_rank(post);
        
        // print result to console
        // printf("Rank of true source: %d (%d)\n", (int) g.rank, g.nsources);
        
        // store rank of true source
        o->qmap[nexp][tinf_idx][i] = g.qmap;
        o->ranks[nexp][tinf_idx][i] = g.rank;
        o->nsources[nexp][tinf_idx][i] = g.nsources;
        
        // compute posterior
        posterior(post);
        
        // compute weighted node predictions
        weighted_node_predictions(tinf_idx, post);
        
        // query next node
        next_node = fn(tinf_idx, post, &evid);
        
        // check if queried node is not NONE
        if (next_node != NONE) {
            // update evidence with new node
            update_evidence(next_node, &evid);
            // set trigger to 1 if evidence is updated
            trigger = 1;
        }
        // else, set trigger to 0 (since we don't need to update likelihood)
        else trigger = 0;
        
    }    
    
    // free memory
    free(post->log_liks); free(post->posterior); free(evid.nodes);
    
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// function to perform querying-inference iteration

void run_querying (unsigned int nexp, POSTERIOR *post, OUTPUTQ *o1, OUTPUTQ *o2, OUTPUTQ *o3, OUTPUTQ *o4, OUTPUTQ *o5) {
    
    // declare integers
    unsigned int i;
    
    // loop over inference times
    for (i = 0; i < g.ntinf; i++) {
        
        // draw first observed node
        draw_first_node(i);
        
        // run inference with different querying strategies
        // as many calls to function as there are querying strategies
        infer_query(nexp, i, post, maxp, o1);
        infer_query(nexp, i, post, ucty, o2);
        infer_query(nexp, i, post, akld, o3);
        infer_query(nexp, i, post, onehop, o4);
        infer_query(nexp, i, post, randq, o5);
        
    }

}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -