#pragma once

#include <list>
#include <vector>
#include <algorithm>
#include <queue>
#include <limits>

#include <iostream>

namespace ln
{
    template<typename value_type>
    class simple_mcf
    {
        static constexpr int alfa = 4;
        int push_count{}, relabel_count{};
        private:
        // graph
        const std::vector< std::pair<int,int> > edge_list;
        const std::vector<int> pos;
        const std::vector< int > adj;
        
        // flow
        std::vector< value_type > residual_capacity;
        
        // cost
        std::vector<value_type> reduced_cost;
        std::vector<value_type> potential,excess;
        
        // path find
        std::vector<int> distance,distance_freq;
        std::vector<int> parent;
        
        // cost-scaling
        value_type epsilon{};
        std::queue<int> active;
        
        static long long int lower_bound_power2(long long int n)
        {
            if(n<=2) return n;
            while(n != (n & -n))
                n -= (n & -n);
            return n;
        }
        
        static auto create_pos(const int N, const std::vector<std::pair<int,int>>& edges)
        {
            std::vector<int> pos(N+1,0);
            for(auto i=0UL;i<edges.size();++i)
            {
                auto [a,b] = edges[i];
                ++pos[a+1];
                ++pos[b+1];
            }
            for(auto i=1UL;i<pos.size();++i)
            {
                pos[i] += pos[i-1];
            }
            return pos;
        }
        
        static auto create_adjacency(const std::vector<int>& pos, const std::vector<std::pair<int,int>>& edges)
        {
            std::vector< int > adj(edges.size()*2);
            std::vector< int > pt(pos.size(),0);
            for(auto i=0UL;i<edges.size();++i)
            {
                auto [a,b] = edges[i];
                adj[ pos[a] + pt[a] ] = 2*i; ++pt[a];
                adj[ pos[b] + pt[b] ] = 2*i+1; ++pt[b];
            }
            return adj;
        }
        void relabel(const int x)
        {
            const value_type relabel_max = std::numeric_limits<value_type>::max();
            value_type delta = relabel_max;
            
            const int xbeg=pos[x],xend=pos[x+1];
            for(int i=xbeg;i<xend;++i)
            {
                int e=adj[i];
                const auto [a,b] = edge_list[e/2];
                const auto next = (e&1) ? a : b;
                
                const auto rcost = reduced_cost[e] + potential[x] - potential[next];
                const auto rcap = residual_capacity[e];
                
                if(rcap>0 && rcost>=0)
                    delta = std::min(delta,rcost+epsilon);
            }
            delta = delta==relabel_max ? epsilon : delta;
            
            ++ relabel_count;
            // reduced_cost[e] = cost[e] + potential[src] - potential[dest]
            
            // std::cerr << "relabel " << x << " value " << delta << '\n';
            potential[x] -= delta;
            // const int xbeg=pos[x],xend=pos[x+1];
            // for(int i=xbeg;i<xend;++i)
            // {
            //     int e=adj[i];
            //     reduced_cost[e] -= delta;
            //     reduced_cost[e^1] += delta;
            // }
        }
        void push_flow(const int arc, value_type delta)
        {
            ++push_count;
            // std::cerr << "push flow in arc " << arc << " value " << delta << '\n';
            residual_capacity[arc]-=delta;
            residual_capacity[arc^1]+=delta;
            
            auto [src,dst] = edge_list[arc/2];
            // if(arc & 1) std::swap(src,dst);
           
            value_type d = (arc & 1) ? 0-delta : delta;
           
            // std::cerr << "from " << src << " to " << dst << '\n';
            excess[src] -= d;
            excess[dst] += d;
        }
            
        void initialize_path_find (const int /*Source*/, const int Dest)
        {
            const int my_INFINITY = std::numeric_limits<int>::max();
            std::fill(parent.begin(),parent.end(),-1);
            std::fill(distance.begin(),distance.end(),my_INFINITY);
            std::fill(distance_freq.begin(),distance_freq.end(),0);
            
            std::queue<int> q;
            distance[Dest]=0;
            
            q.push(Dest);
            
            while(!q.empty())
            {
                auto n = q.front();
                q.pop();
                
                const int xbeg=pos[n],xend=pos[n+1];
                for(int i=xbeg;i<xend;++i)
                {
                    auto e = adj[i];
                    auto arc = e^1;
                    if( residual_capacity[arc]>0 ) 
                    {
                        auto [a,b] = edge_list[e/2];
                        auto next = (e&1) ? a : b;
                        
                        value_type dnew = distance[n] + 1;
                        
                        if(distance[next]==my_INFINITY)
                        {
                            distance[next] = dnew;
                            distance_freq[dnew]++;
                            q.push(next);
                        }
                    }
                }
            }
        }
        
        bool path_find (
            const int Source, const int Dest)
        {
            parent[Dest]=-1;
            const int N = pos.size()-1;
            for(auto current = Source;
                distance[Source]<N && current!=Dest;)
            {
               // advance
               bool found_next=false;
               
               {
                   const int xbeg=pos[current],xend=pos[current+1];
                   for(int i=xbeg;i<xend;++i)
                   {
                        int e = adj[i];
                        auto [a,b] = edge_list[e/2];
                        auto next = (e&1) ? a : b;
                        if(residual_capacity[e]>0 && distance[current]==distance[next]+1)
                        {
                            found_next = true;
                            parent[next] = e;
                            current = next;
                            break;
                        }
                   }
               }
               if(found_next) continue; // advance success
               
               // relabel
               int min_dist = N+10;
               
               {
                   const int xbeg=pos[current],xend=pos[current+1];
                   for(int i=xbeg;i<xend;++i)
                   {
                        int e= adj[i];
                        auto [a,b] = edge_list[e/2];
                        auto next = (e&1) ? a : b;
                        if(residual_capacity[e]>0)
                        {
                            min_dist= std::min(min_dist,distance[next]);
                        }
                   }
               }
               {
                    const int new_dist = min_dist+1;
                    const int old_dist = distance[current];
                    distance[current] = new_dist;
                    if(new_dist<static_cast<int>(distance_freq.size()))
                        distance_freq[new_dist]++;
                    distance_freq[old_dist]--;
                    if(distance_freq[old_dist]==0)
                        break;
               }
               
               // retreat
               if(parent[current]>=0)
               {
                    auto e = parent[current];
                    auto [a,b] = edge_list[e/2];
                    auto prev = (e&1) ? b : a;
                    current = prev;
               }
            }
            return parent[Dest]>=0;
        }
        
        value_type mf_solve(const int Source, const int Dest)
        {
            const value_type my_INFINITY = std::numeric_limits<value_type>::max();
            initialize_path_find(Source,Dest);
            value_type sent=0;
            int cycle=0;
            while(1)
            {
                cycle++;
                // std::cerr << "simple: maxflow cycle " << cycle << '\n';
                bool found = path_find(Source,Dest);
                
                if(!found) break;
                // std::cerr << "path found!\n"; 
                
                value_type k = my_INFINITY;
                for(int current=Dest;parent[current]>=0;)
                {
                    auto e = parent[current];
                    k = std::min(k,residual_capacity[e]);
                    
                    auto [a,b] = edge_list[e/2];
                    auto prev = (e&1) ? b : a;
                    current = prev;
                }
                // std::cerr << " with k = " << k << '\n';
                for(int current=Dest;parent[current]>=0;)
                {
                    auto e = parent[current];
                    
                    // std::cerr << "augment flow in " << e << " value " << k<< '\n';
                   
                    residual_capacity[e] -= k;
                    residual_capacity[e^1] += k;
                    
                    auto [a,b] = edge_list[e/2];
                    auto prev = (e&1) ? b : a;
                    current = prev;
                }
                sent += k;
            }
            return sent;
        }
        
        bool is_admissible(const int arc)const
        {
            const auto [a,b] = edge_list[arc/2];
            const auto next = (arc&1) ? a : b;
            const auto current = (arc&1) ? b : a;
            
            const auto rcost = reduced_cost[arc] + potential[current] - potential[next];
            const auto rcap = residual_capacity[arc];
            
            return (rcost<0 && rcap>0);
        }
        
        bool look_ahead(const int next, const int arc)
        /*
            Inspect the node next if it can receive flow without relabeling
        */
        {
            if(excess[next]<0) return true;
            const int xbeg=pos[next],xend=pos[next+1];
            for(int i=xbeg;i<xend;++i)
            {
                int e=adj[i];
                if(is_admissible(e))
                    return true;
            }
            // node next needs relabeling
            relabel(next);
            return is_admissible(arc);
        }
       
        void examine_node(const int current)
        {
            // std::cerr << "examine node " << current << "\n";
            while(excess[current]>0)
            {
                
                const int xbeg=pos[current],xend=pos[current+1];
                for(int i=xbeg;i<xend;++i)
                {
                    int e=adj[i];
                    const auto [a,b] = edge_list[e/2];
                    const auto next = (e&1) ? a : b;
                    
                    const auto rcost = reduced_cost[e] + potential[current] - potential[next];
                    const auto rcap = residual_capacity[e];
                    
                    // if(rcost<0 && rcost>=-epsilon && rcap>0)
                    if(rcost<0 && rcap>0)
                    {
                        if(!look_ahead(next,e)) continue;
                        
                        // pushed = true;
                        auto d = std::min(excess[current],rcap);
                        
                        push_flow(e,d);
                        
                        if(excess[next]>0)
                        {
                            active.push(next);
                            // std::cerr << "node " << current << " becomes active\n";
                        }
                        if(excess[current]<=0)
                        {
                            // active.erase(current);
                            break;
                            // std::cerr << "node " << current << " no longer active\n";
                        }
                    }
                }
                
                if(excess[current]>0)
                    relabel(current);
            }
        }
                
        bool price_refinement(const value_type eps)
        {
            // std::cerr << "trying price refinement\n";
            // "Price-refinement"
            // check if we can find a price potential such that the current flow is epsilon optimal
            // already, ie. if x[e]>0 then reduced_cost[e]>=-epsilon
            const int N = pos.size()-1;
            const int M = adj.size()/2;
            
            
            // Bellman-Ford
            const value_type INF = std::numeric_limits<value_type>::max();
            std::vector<value_type> dist(N,0);
            
            for(int i=0;i<=N;++i)
            {
                bool updates = false;
                for(int e=0;e<M;++e)
                {
                    const auto [a,b] = edge_list[e];
                   
                    // arc
                    int e1 = e*2;
                    const auto rcost = reduced_cost[e1] + potential[a] - potential[b];
                    if(residual_capacity[e1]>0)
                    {
                        const auto dnew = dist[a]+reduced_cost[e1] + eps;
                        if(dist[b]>dnew)
                        {
                            updates=true;
                            dist[b]=dnew;
                        }
                    }
                    // dual
                    int e2 = e*2+1;
                    if(residual_capacity[e2]>0)
                    {
                        const auto dnew = dist[b]-rcost + eps;
                        if(dist[a]>dnew)
                        {
                            updates = true;
                            dist[a] = dnew;
                        }
                    }
                }
                if(!updates)
                    break;
                
                if(updates && i==N)
                {
                    // std::cerr << "price refinement failed, we hit negative cycle\n";
                    return false;
                }
            }
            
            for(int i=0;i<N;++i)
                potential[i] = dist[i];
            
            // check-condition
            // for(int e=0;e<2*M;++e)
            // if(residual_capacity[e]>0 && reduced_cost[e]+eps<0)
            // {
            //     assert(false);
            //     // return false;
            // }
            
            return true;
        }
        void set_relabel_alternative()
        // Set-relabel according to Bunnagel-Korte-Vygen 1998, 
        // https://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.84.9897
        {
            const int node_size = pos.size()-1;
            const int N = node_size;
            const int INF = std::numeric_limits<int>::max();
            
            value_type remaining_excess=0;
            std::vector< std::list<int> > buckets(1);
            std::vector< int > bucket_value(N,INF);
            for(int node=0;node<N;++node)
            {
                if(excess[node]<0)
                {
                    buckets[0].push_back(node);
                    bucket_value[node]=0;
                }
            }
            
            for(int b=0;b<buckets.size();++b)
            {
                for(auto it = buckets[b].begin();it!=buckets[b].end();++it)
                {
                    const auto node = *it;
                    remaining_excess -= excess[node];
                    
                    const int xbeg=pos[node],xend=pos[node+1];
                    for(int i=xbeg;i<xend;++i)
                    {
                        int e=adj[i];
                        const auto [a,b] = edge_list[e/2];
                        const auto next = (e&1) ? a : b;
                        
                        const auto opposite_arc = e^1;
                        
                        const auto rc = residual_capacity[opposite_arc];
                        const auto rw = 0-(reduced_cost[e]+potential[node]-potential[next]);
                        
                        if(bucket_value[node]<=b || rc<=0)continue;
                        
                        const int d =1+ rw/epsilon;
                        if(bucket_value[node]>d)
                        {
                            bucket_value[node]=d;
                            buckets[d].push_back(node);
                        }
                    }
                }
                if(remaining_excess==0)
                    break;
            }
            
            for(int node=0;node<N;++node)
            {
                const int d = std::min(bucket_value[node],int(buckets.size()-1));
                potential[node] -= epsilon*d;
            }
        }
        void set_relabel()
        // Copyright 2010-2022 Google LLC
        // Set-relaber according to Ortools library, 
        // https://github.com/google/or-tools/blob/stable/ortools/graph/min_cost_flow.cc
        {
            const int node_size = pos.size()-1;
            const int N = node_size;
            std::vector<int> bfs_queue;
            std::vector<bool> node_in_queue(node_size,false);
            
            const value_type kMinCostValue = std::numeric_limits<value_type>::min();
            std::vector<value_type> min_non_admissible_potential(node_size,kMinCostValue);
            std::vector<int> nodes_to_process;
            
            value_type remaining_excess=0;
            
            for(int node=0;node<N;++node)
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
                    
                    const int xbeg=pos[node],xend=pos[node+1];
                    for(int i=xbeg;i<xend;++i)
                    {
                        int e=adj[i];
                        const auto [a,b] = edge_list[e/2];
                        const auto next = (e&1) ? a : b;
                        
                        const auto opposite_arc = e^1;
                        
                        const auto rc = residual_capacity[opposite_arc];
                        const auto rw = 0-(reduced_cost[e]+potential[node]-potential[next]);
                        
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
                                potential[node] - reduced_cost[opposite_arc]
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
                potential_delta = max_potential_diff-epsilon;
                
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
            for(int node=0;node<N;++node)
            {
                if(!node_in_queue[node])
                {
                    potential[node] += potential_delta;
                }
            }
        }
        value_type mcf_solve(const int S, const int T)
        {
            push_count = relabel_count = 0;
            // std::cerr << "simple: maxflow starting\n";
            const value_type maxflow = mf_solve(S,T);
            
            // std::cerr << "simple: maxflow done = " << maxflow << '\n';
            
            epsilon = 0;
            const int N = pos.size()-1;
            const int M = edge_list.size();
            for(auto& w : reduced_cost)
            {
                w *= N;
                epsilon = std::max(epsilon,w);
            }
            // maxC = lower_bound_power2(maxC);
            
            int cycle=0;
            do{
                cycle ++;
                epsilon = std::max(value_type(1),epsilon/alfa);
                // std::cerr << "simple: capacity scaling cycle  " << cycle << '\n';
                // std::cerr << "simple: max-Cost  " << maxC << '\n';
                
                // if(price_refinement(maxC))
                //    continue;
                    
                // std::cerr << "refine\n";
                
                // refine
                for(int e=0;e<2*M;++e)
                {
                    const auto [a,b] = edge_list[e/2];
                    const auto next = (e&1) ? a : b;
                    const auto current = (e&1) ? b : a;
                    
                    const auto rcost = reduced_cost[e] + potential[current] - potential[next];
                    
                    if(rcost<0 && residual_capacity[e]>0)
                    {
                        push_flow(e,residual_capacity[e]);
                    }
                    // this also does the trick of making the flow = 0 for arcs with
                    // reduced_weight>0, 
                }
                for(int n=0;n<N;++n)
                if(excess[n]>0)
                {
                    active.push(n);
                }
                while(!active.empty())
                {
                    if(relabel_count>=N)
                    {
                        relabel_count=0;
                        set_relabel();
                    }
                    auto current = active.front();
                    active.pop();
                    
                    examine_node(current);
                    
                }
            }while(epsilon>1);
            return maxflow;
        }
        
        public:
        simple_mcf(
            const int N_nodes,
            const std::vector<std::pair<int,int>>& edges,
            const std::vector<value_type>& cap,
            const std::vector<value_type>& cost):
                edge_list{edges}, 
                pos{create_pos(N_nodes,edges)},
                adj{create_adjacency(pos,edges)},
                residual_capacity(edges.size()*2),reduced_cost(edges.size()*2),
                potential(N_nodes,0),
                excess(N_nodes,0),
                distance(N_nodes,0),
                distance_freq(N_nodes+1,0),
                parent(N_nodes,-1)
        {
            const int M = edges.size();
            for(int i=0;i<M;++i)
            {
                residual_capacity[2*i] = cap[i];
                residual_capacity[2*i+1]=0;
                
                reduced_cost[2*i] = cost[i];
                reduced_cost[2*i+1] = 0-cost[i];
            }
        }
        
        auto solve(const int S, const int T)
        {
            const int M = edge_list.size();
            mcf_solve(S,T);
            
            std::vector<value_type> flow(edge_list.size());
            for(int e=0;e<M;++e)
            {
                const auto arc = e*2;
                const auto dual = arc+1;
                flow[e] = residual_capacity[dual];
            }
            
            return flow;
        }
    };
}
