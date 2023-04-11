#include <mincostflow/mincostflow.hpp>
#include <iostream>
#include <chrono>

#ifdef USE_ORTOOLS
#include <ortools/graph/min_cost_flow.h>
#include <ortools/graph/max_flow.h>
#endif

#include <mincostflow/simple.hpp>

typedef long long value_type;
typedef int nodeID_type;
typedef int arcID_type;
typedef ln::digraph<nodeID_type,arcID_type> graph_type;

template<typename solver_t>
struct my_mcf
{
    const std::size_t M;
    graph_type g;
    std::vector<value_type> capacity;
    std::vector<value_type> weight;
    
    my_mcf(
        const std::size_t N_nodes,
        const std::vector<std::pair<nodeID_type,nodeID_type>>& edges,
        const std::vector<value_type>& cap,
        const std::vector<value_type>& cost):
            M{edges.size()}
    {
        for(nodeID_type i=0;i<static_cast<nodeID_type>(N_nodes);++i)
            g.add_node(i);
        
        
        for(arcID_type e=0;e<static_cast<arcID_type>(M);++e)
        {
            g.add_arc(edges[e].first,edges[e].second,e);
        }
        
        set_capacity(cap);
        set_cost(cost);
    }
    
    auto solve(nodeID_type S, nodeID_type T)
    {
        std::vector<value_type> flow(M);
        
        solver_t mcf;
        mcf.solve(g,g.get_node(S),g.get_node(T),weight,capacity);
        
        for(arcID_type e=0;e<static_cast<arcID_type>(M);++e)
        {
            auto arc = g.get_arc(e);
            auto dual = g.arc_dual(arc);
            flow[e] = capacity[dual];
        }
        
        return flow;
    }
    
    void set_capacity(const std::vector<value_type>& cap)
    {
        capacity.resize(g.max_num_arcs());
        for(arcID_type i=0;i<static_cast<arcID_type>(cap.size());++i)
        {
            auto arc = g.get_arc(i);
            auto dual = g.arc_dual(arc);
            
            capacity[arc] = cap[i];
            capacity[dual]=0;
        }
    }
    void set_cost(const std::vector<value_type>& cost)
    {
        weight.resize(g.max_num_arcs());
        for(arcID_type i=0;i<static_cast<arcID_type>(cost.size());++i)
        {
            auto arc = g.get_arc(i);
            auto dual = g.arc_dual(arc);
            
            weight[arc] = cost[i];
            weight[dual]= -cost[i];
        }
    }
};

template<typename solver_type>
std::pair<value_type,value_type> solve_generic(
    const std::size_t N_nodes,
    nodeID_type S, nodeID_type T,
    const std::vector<std::pair<nodeID_type,nodeID_type>>& edges,
    const std::vector<value_type>& cap,
    const std::vector<value_type>& cost,
    const std::string tname)
{
    const std::size_t M = edges.size();
    assert(cap.size()==M);
    assert(cost.size()==M);
    
    // initialize the solver 
    solver_type solver(N_nodes,edges,cap,cost);
   
   
   // solve
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<value_type> flow = solver.solve(S,T);
   
    auto stop = std::chrono::high_resolution_clock::now();
    
    std::cout << tname << " " <<
    std::chrono::duration_cast<std::chrono::microseconds>(stop-start) .count()
    << std::endl;
    
    // check the capacity constraints
    for(arcID_type e=0;e<static_cast<arcID_type>(M);++e)
    {
        assert(flow[e]>=0 && flow[e]<=cap[e]);
    }
    
    // get the balance
    std::vector<value_type> balance(N_nodes,0);
    for(arcID_type e=0;e<static_cast<arcID_type>(M);++e)
    {
        auto [a,b] = edges[e];
        auto f = flow[e];
        balance[a] -= f;
        balance[b] += f;
    }
    
    // check constraints on the balance
    assert(balance[S]==-balance[T]);
    assert(balance[T]>=0);
    for(nodeID_type i=0;i<static_cast<nodeID_type>(N_nodes);++i)
        if(i!=S && i!=T)
            assert(balance[i]==0);
    
    value_type mcost = 0;
    value_type mflow = balance[T];
    
    // compute the cost
    for(arcID_type e=0;e<static_cast<arcID_type>(M);++e)
    {
        mcost += flow[e] * cost[e];
    }
    
    return {mflow,mcost};
}

#ifdef USE_ORTOOLS
struct ortool_mcf
{
    const std::size_t M;
    operations_research::SimpleMaxFlow max_flow;
    operations_research::SimpleMinCostFlow mincost_flow;
    
    ortool_mcf(
        const std::size_t N_nodes,
        const std::vector<std::pair<nodeID_type,nodeID_type>>& edges,
        const std::vector<value_type>& cap,
        const std::vector<value_type>& cost):
            M{edges.size()}
    {
        for(arcID_type e=0;e<static_cast<arcID_type>(M);++e)
        {
            auto [a,b] = edges[e];
            max_flow.AddArcWithCapacity(a,b,cap[e]);
            mincost_flow.AddArcWithCapacityAndUnitCost(a,b,cap[e],cost[e]);
        }
        for(nodeID_type i=0;i<static_cast<nodeID_type>(N_nodes);++i)
        {
            mincost_flow.SetNodeSupply(i,0);
        }
    }
    
    auto solve(nodeID_type S, nodeID_type T)
    {
        /* int status = */ max_flow.Solve(S,T);
        const value_type Flow = max_flow.OptimalFlow();
        
        mincost_flow.SetNodeSupply(S,Flow);
        mincost_flow.SetNodeSupply(T,-Flow);
        
        /* int min_status = */ mincost_flow.Solve();
        // const value_type Cost = mincost_flow.OptimalCost();
        
        std::vector<value_type> flow(M);
        
        for(arcID_type e=0;e<static_cast<arcID_type>(M);++e)
        {
            flow[e] = mincost_flow.Flow(e);
        }
        
        return flow;
    }
    
};
#endif

int main()
{   
    int N,M,S,T;
    std::cin >> N >> M >> S >> T;
    
    std::vector<std::pair<int,int>> ed_list; 
    std::vector<value_type> capacity;
    std::vector<value_type> weight;
    
    for(int e=0;e<M;++e)
    {
        int a,b,wei,cap;
        std::cin>>a>>b>>cap>>wei;
        
        ed_list.push_back({a,b});
        capacity.push_back(cap);
        weight.push_back(wei);
    }
    
    
    auto [flow_0,cost_0] =
    solve_generic<ln::simple_mcf<value_type>>(N,S,T,ed_list,capacity,weight,"Cost-scaling-simple");
    
    {
        auto [flow,cost] = solve_generic<
                my_mcf<  ln::mincostflow_EdmondsKarp<value_type,ln::shortestPath_FIFO<value_type>> >
                 >(N,S,T,ed_list,capacity,weight,"Augmenting-path");
        assert(flow_0==flow && cost_0==cost);
    }
    {
        auto [flow,cost] = 
                solve_generic<
                    my_mcf< 
                        ln::mincostflow_PrimalDual<
                            ln::shortestPath_Dijkstra<value_type>,
                            ln::maxflow_augmenting_path<
                                value_type,
                                ln::pathSearch_labeling
                            >
                        >
                    >
                >(N,S,T,ed_list,capacity,weight,"Primal-dual");
        assert(flow_0==flow && cost_0==cost);
    }
    { 
        auto [flow,cost] = 
                solve_generic<
                    my_mcf< 
                        ln::mincostflow_costScaling<
                            ln::maxflow_augmenting_path<
                                value_type,
                                ln::pathSearch_labeling
                            >
                        >
                    >
                >(N,S,T,ed_list,capacity,weight,"Cost-scaling");
        assert(flow_0==flow && cost_0==cost);
    }
    #ifdef USE_ORTOOLS
    {
        auto [flow,cost] = solve_generic<ortool_mcf>(N,S,T,ed_list,capacity,weight,"Ortools");
        assert(flow_0==flow && cost_0==cost);
    }
    #endif
    return 0;
}

