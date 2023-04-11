#pragma once

#include <mincostflow/maxflow.hpp>
#include <mincostflow/scope_guard.hpp>

#include <queue>

namespace ln
{
    template<typename T, typename path_optimizer_type>
    class mincostflow_EdmondsKarp : public maxflow_base<T>
    {
        public:
        using base_type = maxflow_base<T>;
        using value_type = typename base_type::value_type;    
        using node_pos_t = typename base_type::node_pos_t;
        using arc_pos_t = typename base_type::arc_pos_t;
        using base_type::flow_at;
        using base_type::INFINITY;
        
        template<typename graph_t>
        value_type solve(
            const graph_t& g,
            const node_pos_t Source, const node_pos_t Dest,
            const std::vector<value_type>& weight,
                  std::vector<value_type>& residual_cap
            )
        // augmenting path
        {   
            value_type sent =0 ;
            path_optimizer_type path_opt;
            
            while(true)
            {
                path_opt.solve(
                    g,
                    Source,
                    weight,
                    // edge is valid if
                    [&residual_cap](arc_pos_t e){
                        return residual_cap.at(e)>0;
                    });
                
                if(! path_opt.is_reacheable(Dest))
                    break;
                
                auto path = path_opt.get_path(g,Dest);
                
                value_type k = INFINITY;
                for(auto e : path)
                {
                    k = std::min(k,residual_cap.at(e));
                }
                
                for(auto e: path)
                {
                    residual_cap[e] -= k;
                    residual_cap[g.arc_dual(e)] += k;
                } 
                
                sent += k;
            }
            return sent;
        }
        
        mincostflow_EdmondsKarp()
        {}
    };
    
    
    template<typename path_optimizer_type, typename maxflow_type>
    class mincostflow_PrimalDual : public maxflow_type
    {
        public:
        using base_type = maxflow_type;
        using value_type = typename base_type::value_type;    
        using node_pos_t = typename base_type::node_pos_t;
        using arc_pos_t = typename base_type::arc_pos_t;
        using base_type::flow_at;
        using base_type::INFINITY;
        
        
        template<typename graph_t>
        value_type solve(
            const graph_t& g,
            const node_pos_t Source, const node_pos_t Dest,
            const std::vector<value_type>& weight,
                  std::vector<value_type>& residual_cap
            )
        {   
            std::vector<value_type> reduced_weight = weight;
            
            value_type sent =0 ;
            path_optimizer_type path_opt;
            
            while(true)
            {
                path_opt.solve(
                    g,
                    Source,
                    reduced_weight,
                    // edge is valid if
                    [&residual_cap](arc_pos_t e) -> bool
                    {
                        return residual_cap.at(e)>0;
                    });
                    
                if(! path_opt.is_reacheable(Dest))
                    break;
                    
                const auto& distance{path_opt.distance};
                
                for(auto e : g.arcs())
                {
                
                    auto [a,b] = g.arc_ends(e);
                    if(distance[a]<INFINITY && distance[b]<INFINITY)
                    {
                        reduced_weight[e]       += distance[a]-distance[b];
                    }
                }
                
                
                auto F = base_type::solve(
                    g,
                    Source,Dest,
                    residual_cap,
                    // admissibility
                    [&reduced_weight](arc_pos_t e)->bool
                    {
                        return reduced_weight[e]==0;
                    });
                
                sent += F;
            }
            return sent;
        }
        
        mincostflow_PrimalDual()
        {}
    };
    
    template<typename path_optimizer_type, typename maxflow_type>
    class mincostflow_capacityScaling : public maxflow_type
    {
        public:
        using base_type = maxflow_type;
        using value_type = typename base_type::value_type;
        using node_pos_t = typename base_type::node_pos_t;
        using arc_pos_t  = typename base_type::arc_pos_t;
        using base_type::flow_at;
        using base_type::INFINITY;
        
        template<typename graph_t>
        value_type solve(
            graph_t& g,
            const node_pos_t Source, const node_pos_t Dest,
            const std::vector<value_type>& weight,
                  std::vector<value_type>& residual_cap)
        {
        
            std::vector<value_type> reduced_weight = weight;
            value_type maxflow{0};
            
            // find the max-flow-anycost
            maxflow = maxflow_type::solve(
                g,Source,Dest,
                residual_cap,
                [](arc_pos_t)->bool{return true;});
            
            value_type cap_flow = lower_bound_power2(maxflow);
            
            std::vector<value_type> excess(g.max_num_nodes(),0);
            
            std::vector<value_type> weight_ex = weight;
            
            auto update_reduced_costs = 
                [&](const std::vector<value_type>& potential)
            {
                for(auto e : g.arcs())
                {
                    auto [src,dst] = g.arc_ends(e);
                    auto p_src = potential.at(src), p_dst = potential.at(dst);
                    
                    p_src = p_src == INFINITY ? 0 : p_src;
                    p_dst = p_dst == INFINITY ? 0 : p_dst;
                    
                    weight_ex.at(e) +=  p_src - p_dst;
                }
            };
            
            auto push_flow = 
                [&](arc_pos_t e,value_type delta)
            {
                    // std::cerr << " push flow at " << e << " delta = " << delta << "\n";
                    auto [src,dst] = g.arc_ends(e);
                    
                    residual_cap[e]-=delta;
                    residual_cap[g.arc_dual(e)]+=delta;
                    
                    excess.at(src) -= delta;
                    excess.at(dst) += delta;
                    
                    // std::cerr << "push " << delta << " over " << e << '\n';
            };
            
            // auto report = 
            // [&]()
            // {
            //     std::cerr << "residual cap + mod. costs\n";
            //     for(auto e : g.arcs())
            //     {
            //         std::cerr << " " << e << " -> " << residual_cap[e] << " " << weight_ex[e] << "\n";
            //     }
            //     std::cerr << "potential + excess\n";
            //     for(auto v : g.nodes())
            //     {
            //         std::cerr << " " << v << " -> " << excess[v] << "\n";
            //     }
            // };
            // 
            // std::cerr << " maxflow = " << maxflow << "\n";
            
            // int cycle=0;
            for(;cap_flow>0;cap_flow/=2)
            {
                // cycle++;
                // std::cerr << "cycle " << cycle << " cap_flow = " << cap_flow << '\n';
                // report();
                
                // saturate edges with negative cost
                for(auto e : g.arcs()) 
                while(residual_cap.at(e)>=cap_flow && weight_ex.at(e)<0)
                {
                    push_flow(e,cap_flow);
                }
                
                path_optimizer_type path_opt;
                
                // build S and T
                std::set<node_pos_t> Sset,Tset;
                for(auto v : g.nodes())
                {
                    if(excess.at(v)>=cap_flow)
                        Sset.insert(v);
                    if(excess.at(v)<=-cap_flow)
                        Tset.insert(v);
                }
                
                const auto multi_source_node = g.new_node();
                excess.resize(g.max_num_nodes());
                excess.at(multi_source_node) = 0;
                const Scope_guard rm_node = [&](){ g.erase_node(multi_source_node);};
                
                
                
                for(auto v : Sset)
                {
                    auto arc1 = g.new_arc(multi_source_node,v);
                    auto arc2 = g.new_arc(v,multi_source_node);
                    
                    g.set_dual(arc1,arc2);
                    
                    weight_ex.resize(g.max_num_arcs());
                    residual_cap.resize(g.max_num_arcs());
                    
                    weight_ex.at(arc1) = 0;
                    residual_cap.at(arc1) = excess.at(v);
                    
                    weight_ex.at(arc2) = 0;
                    residual_cap.at(arc2) = 0;
                    
                    excess.at(multi_source_node) += excess.at(v);
                    excess.at(v) = 0;
                }
                
                const Scope_guard restore_excess = [&]()
                {
                    for(auto e : g.out_arcs(multi_source_node))
                    {
                        auto [src,dst] = g.arc_ends(e);
                        excess.at(dst) = residual_cap.at(e);
                    }
                };
                
                while(!Sset.empty() && !Tset.empty())
                { 
                    path_opt.solve(
                        g, multi_source_node,
                        weight_ex,
                        [cap_flow,&residual_cap](arc_pos_t e)->bool
                        {
                            return residual_cap.at(e)>=cap_flow;
                        }
                    );
                    
                    const auto& distance{path_opt.distance};
                    
                    auto it = std::find_if(Tset.begin(),Tset.end(),
                        [&](node_pos_t v)->bool {
                            return distance.at(v)<INFINITY;
                        });
                    
                    if(it==Tset.end())
                        break;
                    
                    auto dst = *it;
                    
                    
                    // std::cerr << " vertex distance to pivot\n";
                    // for(int v=0;v<Graph.n_vertex();++v)
                    // {
                    //     std::cerr << " " <<v <<" -> " << distance[v]<<"\n";
                    // }
                    
                    update_reduced_costs(distance);
                    
                    auto path = path_opt.get_path(g,dst);
                    for(auto e: path)
                    {
                        // auto [src,dst] = g.arc_ends(e);
                        push_flow(e,cap_flow);
                    }
                    
                    if(excess.at(dst)>-cap_flow)
                        Tset.erase(dst);
                }
            }
            
            return maxflow;
        }
        
        public:
        mincostflow_capacityScaling()
        {}
    };
    
    template<typename maxflow_type>
    class mincostflow_costScaling : public maxflow_type
    {
        public:
        using base_type = maxflow_type;
        using value_type = typename base_type::value_type;
        using node_pos_t = typename base_type::node_pos_t;
        using arc_pos_t  = typename base_type::arc_pos_t;
        using base_type::flow_at;
        using base_type::INFINITY;
        
        private:
        
        static constexpr int alpha = 4; // the factor that divides epsilon at each refinement step
        unsigned int push_count{}, relabel_count{};
        value_type eps{};
            
        std::vector<value_type> potential;
        std::vector<value_type> excess;
        std::vector<value_type> residual_cap;
        std::vector<value_type> scaled_cost;
        std::queue<node_pos_t> active;
        
        template<typename graph_t>
        void relabel(const graph_t& g, const node_pos_t x)
        {
            const value_type delta_max = std::numeric_limits<value_type>::max(); 
            value_type delta = delta_max;
            
            for(auto e : g.out_arcs(x))
            {
                const auto [a,b] = g.arc_ends(e);
                const auto rw = scaled_cost[e] + potential[a] - potential[b];
                const auto rc = residual_cap[e];
                
                if(rc>0 && rw>=0)
                    delta = std::min(delta,rw+eps);
            }
            
            delta = delta==delta_max ? eps : delta ;
            
            // std::cerr << "relabel " << x << " with value " << delta << "\n";
            potential[x] -= delta;
            ++relabel_count;
        }
        
        template<typename graph_t>
        void push(const graph_t& g, const arc_pos_t arc, const value_type delta)
        {
            residual_cap[arc]-=delta;
            residual_cap[g.arc_dual(arc)]+=delta;
                    
            const auto [src,dst] = g.arc_ends(arc);
            // std::cerr << "pushing "<<delta<<" on arc "<<arc<<" ("<<src<<","<<dst<<")\n";
            excess[src] -= delta;
            excess[dst] += delta;
            ++push_count;
        }
        
        template<typename graph_t>
        bool is_admissible(const graph_t& g, const arc_pos_t arc)const
        {
            const auto [a,b] = g.arc_ends(arc);
            const auto rw = scaled_cost[arc] + potential[a] - potential[b];
            const auto rc = residual_cap[arc];
            
            return rw<0 && rc>0;
        }
        
        template<typename graph_t>    
        bool look_ahead(const graph_t& g, const node_pos_t next, const arc_pos_t arc)
        /*
            Inspect the node next if it can receive flow without relabeling
        */
        {
            if(excess[next]<0) return true;
            for(auto e : g.out_arcs(next))
            {
                if(is_admissible(g,e))
                    return true;
            }
            // node next needs relabeling
            relabel(g,next);
            return is_admissible(g,arc);
        }
            
        template<typename graph_t>
        void examine_node(const graph_t& g, const node_pos_t current)
        {
            // std::cerr << "examine node " << current << '\n';
            int cycle = 0;
            while(excess[current]>0)
            {
                cycle++;
                // std::cerr << "examine node " << current << " cycle " << cycle << '\n'; 
                //if(cycle>10)
                //    throw std::runtime_error("ops");
                
                for(auto e : g.out_arcs(current))
                {
                    const auto [a,b] = g.arc_ends(e);
                    const auto rw = scaled_cost[e] + potential[a] - potential[b];
                    const auto rc = residual_cap[e];
                    
                    if(rw<0 && rc>0)
                    {
                        // std::cerr << "admissible arc " << e << '\n';
                        
                        if(!look_ahead(g,b,e)) continue;
                        const auto d = std::min(excess[a],rc);
                        
                        push(g,e,d);
                        
                        if(excess[b]>0)
                        {
                            active.push(b);
                        }
                        if(excess[current]<=0)
                        {
                            break;
                        }
                    }
                }
                
                if(excess[current]>0)
                    relabel(g,current);
            }
        }
        
        template<typename graph_t>
        void set_relabel(const graph_t& g)
        // Copyright 2010-2022 Google LLC
        // Set-relaber according to Ortools library, 
        // https://github.com/google/or-tools/blob/stable/ortools/graph/min_cost_flow.cc
        {
            const int node_size = g.max_num_nodes();
            std::vector<node_pos_t> bfs_queue;
            std::vector<bool> node_in_queue(node_size,false);
            
            const value_type kMinCostValue = std::numeric_limits<value_type>::min();
            std::vector<value_type> min_non_admissible_potential(node_size,kMinCostValue);
            std::vector<node_pos_t> nodes_to_process;
            
            value_type remaining_excess=0;
            
            for(const auto node : g.nodes())
            {
                if(excess[node]<0)
                {
                    bfs_queue.push_back(node);
                    node_in_queue[node]=true;
                    remaining_excess -= excess[node];
                }
            }
            
            value_type potential_delta =0 ;
            unsigned int queue_index=0;
            while(remaining_excess>0)
            {
                for(;queue_index<bfs_queue.size();++queue_index)
                {
                    const auto node = bfs_queue[queue_index];
                    for(const auto e : g.out_arcs(node))
                    {
                        const auto [a,next] = g.arc_ends(e);
                        const auto opposite_arc = g.arc_dual(e);
                        
                        const auto rc = residual_cap[opposite_arc];
                        const auto rw = 0-(scaled_cost[e]+potential[node]-potential[next]);
                        
                        if(node_in_queue[next] || rc<=0)continue;
                        
                        potential[next] += potential_delta;
                        if(rw<0)
                        {
                            remaining_excess -= excess[next];
                            bfs_queue.push_back(next);
                            node_in_queue[next]=true;
                            
                            if(remaining_excess==0)
                            {
                                potential[next] -= potential_delta;
                                break;
                            }
                        }
                        else
                        {
                            potential[next] -= potential_delta;
                            if(min_non_admissible_potential[next]==kMinCostValue)
                            {
                                nodes_to_process.push_back(next);
                            }
                            min_non_admissible_potential[next] = std::max(
                                min_non_admissible_potential[next],
                                potential[node] - scaled_cost[opposite_arc]
                            );
                        }
                    }
                    if(remaining_excess==0)break;
                }
                
                if(remaining_excess==0)break;
                
                value_type max_potential_diff = kMinCostValue;
                for(const auto node : nodes_to_process)
                {
                    if(node_in_queue[node])continue;
                    max_potential_diff =
                        std::max(max_potential_diff,min_non_admissible_potential[node]-potential[node]);
                    if(max_potential_diff==potential_delta)break;
                }
                potential_delta = max_potential_diff-eps;
                
                int index=0;
                for(unsigned int i=0;i<nodes_to_process.size();++i)
                {
                    const auto node = nodes_to_process[i];
                    if(node_in_queue[node])continue;
                    if(potential[node]+potential_delta<min_non_admissible_potential[node])
                    {
                        potential[node]+=potential_delta;
                        bfs_queue.push_back(node);
                        node_in_queue[node]=true;
                        remaining_excess -= excess[node];
                        continue;
                    }
                    nodes_to_process[index]=node;
                    index++;
                }
                nodes_to_process.resize(index);
            }
            
            if(potential_delta==0)return;
            for(const auto node : g.nodes())
            {
                if(!node_in_queue[node])
                {
                    potential[node] += potential_delta;
                }
            }
        }
        
        public:
        
        template<typename graph_t>
        value_type solve(
            const graph_t& g,
            const node_pos_t Source, const node_pos_t Dest,
            const std::vector<value_type>& weight,
                  std::vector<value_type>& input_capacity)
        {
            residual_cap = input_capacity;
            scaled_cost = weight;
            push_count=0;
            relabel_count=0;
            
            value_type maxflow{0};
            
            // find the max-flow-anycost
            maxflow = maxflow_type::solve(
                g,Source,Dest,
                residual_cap,
                [](arc_pos_t)->bool{return true;});
           
            potential.assign(g.max_num_nodes(),0);
            excess.assign(g.max_num_nodes(),0);
            
            
            
            eps = 0;
            const int N = g.num_nodes();
            for(auto e : g.arcs())
            {
                scaled_cost[e] *= N;
                eps = std::max(eps,scaled_cost[e]);
            }
            
            
            int cycle=0;
            do{
                
                eps = std::max(value_type(1),eps/alpha);
                
                cycle++;
                //std::cerr << "cycle " << cycle << " eps = " << eps << '\n';
            
                // improve
                for(auto e : g.arcs())
                {
                    const auto [a,b] = g.arc_ends(e);
                    const auto rw = scaled_cost[e] + potential[a] - potential[b];
                    
                    if(rw<0 && residual_cap[e]>0)
                    {
                        push(g,e,residual_cap[e]);
                    }
                    // this also does the trick of making the flow = 0 for arcs with
                    // scaled_cost>0, 
                }
                for(auto n : g.nodes())
                if(excess[n]>0)
                {
                    active.push(n);
                    // std::cerr << n << " is added to active\n";
                }
                
                int ct=0;
                while(!active.empty())
                {
                    if(relabel_count>=g.num_nodes())
                    {
                        relabel_count=0;
                        set_relabel(g);
                    }
                    
                    ct++;
                    //if(ct>20)throw std::runtime_error("shit");
                
                    auto current = active.front();
                    active.pop();
                    examine_node(g,current);
                }
                 
            }while(eps>1);
            // std::cerr << "done\n";
            input_capacity = residual_cap;
            return maxflow;
        }
        
        public:
        mincostflow_costScaling()
        {}
    };

}
