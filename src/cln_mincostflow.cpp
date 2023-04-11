#include <mincostflow/ln_mincostflow.h>
#include <mincostflow/error.hpp>

extern "C" {

/* Initialize the Network */

MCFNetwork* MCFNetwork_new()
{
    return new MCFNetwork();
}

int MCFNetwork_SetNodeSupply(
    MCFNetwork* g,
    const nodeID_type id, 
    const costflow_value_type supply)
{
    try{
        g->SetNodeSupply(id,supply);
    }catch(ln::mcf_exception& e)
    {
        return e.error_code();
    }catch(...)
    {
        return ln::error_codes::UNKNOWN_ERROR;
    }
    return ln::error_codes::NO_ERROR;
}

int MCFNetwork_AddArc(
    MCFNetwork* g,
    const nodeID_type tail, 
    const nodeID_type head, 
    const arcID_type id,
    const costflow_value_type capacity,
    const costflow_value_type cost,
    std::uint64_t * p_arc_index)
{
    try{
        *p_arc_index = g->AddArc(tail,head,id,capacity,cost);
    }catch(ln::mcf_exception & e)
    {
        return e.error_code();
    }catch(...)
    {
        return ln::error_codes::UNKNOWN_ERROR;
    }
    return ln::error_codes::NO_ERROR;
}

/* Compute */

int MCFNetwork_Solve(MCFNetwork* g)
{
    try
    {
        g->Solve();
    }catch(ln::mcf_exception & e)
    {
        return e.error_code();
    }
    catch(...)
    {
        return ln::error_codes::UNKNOWN_ERROR;
    }
    return ln::error_codes::NO_ERROR;
}

int MCFNetwork_Solve_by_AugmentingPaths(MCFNetwork* g)
{
    try
    {
        g->Solve_by_AugmentingPaths();
    }catch(ln::mcf_exception & e)
    {
        return e.error_code();
    }
    catch(...)
    {
        return ln::error_codes::UNKNOWN_ERROR;
    }
    return ln::error_codes::NO_ERROR;
}

int MCFNetwork_Solve_by_CostScaling(MCFNetwork* g)
{
    try
    {
        g->Solve_by_CostScaling();
    }catch(ln::mcf_exception& e)
    {
        return e.error_code();
    }
    catch(...)
    {
        return ln::error_codes::UNKNOWN_ERROR;
    }
    return ln::error_codes::NO_ERROR;
}

/* Extract information */         

int MCFNetwork_Flow(
    const MCFNetwork* g,
    const arcID_type a,
    costflow_value_type * p_value)
{
    try{
        *p_value = g->Flow(a);
    }catch(ln::mcf_exception& e)
    {
        return e.error_code();
    }catch(...)
    {
        return ln::error_codes::UNKNOWN_ERROR;
    }
    return ln::error_codes::NO_ERROR;
}

int MCFNetwork_CountFlow(
    MCFNetwork const * g,
    std::uint64_t * p_value)
{
    try{
        *p_value = g->CountFlow();
    }catch(ln::mcf_exception& e)
    {
        return e.error_code();
    }catch(...)
    {
        return ln::error_codes::UNKNOWN_ERROR;
    }
    return ln::error_codes::NO_ERROR;
}

int MCFNetwork_GetFlow(
    MCFNetwork const * g,
    uint64_t * const index_list, 
    costflow_value_type * const flow_list)
{
    try{
        g->GetFlow(index_list,flow_list);
    }catch(ln::mcf_exception& e)
    {
        return e.error_code();
    }catch(...)
    {
        return ln::error_codes::UNKNOWN_ERROR;
    }
    return ln::error_codes::NO_ERROR;
}


/* Update the Network */

int MCFNetwork_UpdateArc(
    MCFNetwork* g,
    const arcID_type id,
    const costflow_value_type capacity,
    const costflow_value_type cost)
{
    try{
        g->UpdateArc(id,capacity,cost);
    }catch(ln::mcf_exception& e)
    {
        return e.error_code();
    }catch(...)
    {
        return ln::error_codes::UNKNOWN_ERROR;
    }
    return ln::error_codes::NO_ERROR;
}

int MCFNetwork_RemoveArc(
    MCFNetwork*g, 
    const arcID_type id)
{
    try{
        g->RemoveArc(id);
    }catch(ln::mcf_exception& e)
    {
        return e.error_code();
    }catch(...)
    {
        return ln::error_codes::UNKNOWN_ERROR;
    }
    return ln::error_codes::NO_ERROR;
}

int MCFNetwork_forget(MCFNetwork *g)
{
    try{
        g->forget();
    }catch(ln::mcf_exception& e)
    {
        return e.error_code();
    }catch(...)
    {
        return ln::error_codes::UNKNOWN_ERROR;
    }
    return ln::error_codes::NO_ERROR;
}

/* Cleanup */    

void MCFNetwork_free(MCFNetwork* g)
{
    delete g;
}

}
