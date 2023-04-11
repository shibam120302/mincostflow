#pragma once

#ifdef __cplusplus

#include <mincostflow/ln_mincostflow.hpp>

extern "C" {

using ln::costflow_value_type;

using ln::nodeID_type;
using ln::arcID_type;
using ln::graph_type;
using ln::MCFNetwork;
    
    
#else

#include <stdint.h>

typedef int64_t costflow_value_type;

struct nodeID_type;
struct arcID_type;
struct graph_type;
struct MCFNetwork;

#endif

/* Initialize the Network */

MCFNetwork* MCFNetwork_new();

int MCFNetwork_SetNodeSupply(
    MCFNetwork* g,
    const nodeID_type id, 
    const costflow_value_type supply);

int MCFNetwork_AddArc(
    MCFNetwork* g,
    const nodeID_type tail, 
    const nodeID_type head, 
    const arcID_type id,
    const costflow_value_type capacity,
    const costflow_value_type cost,
    uint64_t * p_arc_index);

/* Compute */

int MCFNetwork_Solve(MCFNetwork* g);

// These two API entries are for benchmarking only
int MCFNetwork_Solve_by_AugmentingPaths(MCFNetwork* g);
int MCFNetwork_Solve_by_CostScaling(MCFNetwork* g);


/* Extract information */         

int MCFNetwork_Flow(
    MCFNetwork const * g,
    const arcID_type a,
    costflow_value_type* p_value);

int MCFNetwork_CountFlow(
    MCFNetwork const* g,
    uint64_t* p_n);

int MCFNetwork_GetFlow(
    MCFNetwork const * g,
    uint64_t * const index_list, 
    costflow_value_type * const flow_list);

/* Update the Network */

int MCFNetwork_UpdateArc(
    MCFNetwork* g,
    const arcID_type id,
    const costflow_value_type capacity,
    const costflow_value_type cost);
            
int MCFNetwork_RemoveArc(
    MCFNetwork*g, 
    const arcID_type id);

int MCFNetwork_forget(MCFNetwork *g);
    
/* Cleanup */    
    
void MCFNetwork_free(MCFNetwork* g);


#ifdef __cplusplus
}
#endif
