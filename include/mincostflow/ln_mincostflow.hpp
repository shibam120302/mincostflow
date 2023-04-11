#pragma once

#include <mincostflow/mincostflow.hpp>
#include <mincostflow/hash_combine.hpp>
#include <mincostflow/simple.hpp>
#include <cstdint>
#include <vector>
#include <array>

#include <mincostflow/error.hpp>

namespace ln
{
    struct nodeID_type
    {
        // Node id is the public key
        
        static constexpr int KEY_SIZE = 33;
        typedef std::uint8_t value_type;
        
        value_type k[KEY_SIZE];
        
        bool operator == (const nodeID_type& that)const
        {
            bool result=true;
            for(auto i=0UL;i<KEY_SIZE;++i)
                result = result && k[i]==that.k[i];
            return result;
        }
    };
    
    
    struct arcID_type
    {
        // arc ID is the combination of the short_channel_id and an index called `part`.
        // `part` can be anything that breaks the degeneracy of the map channel to arcs.
        // there is a degeneracy of channels to arcs because channels come with two possible
        // directions, and for each direction we have an arc with convex cost wich Pickhardtpayments
        // decomposes into several linear-cost arcs (this is the piecewise linear approximation of
        // the cost).
        // One possible degeneracy breaking scheme could be to enumerate each arc piece with a
        // number `index` and define another single bit number `b` to be 0 if the end nodes public keys
        // satisfy A<B or `b`=1 if A>B. Then `part` = `index` << 1  | b.
        
        std::uint64_t short_channel_id;
        std::uint8_t part; // 1 bit for the direction 7 bits for the linearization part
        
        bool operator == (const arcID_type& that)const
        {
            return short_channel_id==that.short_channel_id && part==that.part;
        }
    };
}

namespace std
{
    template<>
    struct hash<ln::nodeID_type>
    {
        // hash function needed to store nodeID_type into std::unordered_map
        
        std::size_t operator()(const ln::nodeID_type& id)const
        {
            std::hash<ln::nodeID_type::value_type> hasher;
            std::size_t seed{0};
            for(auto x : id.k)
                seed = ln::hash_combine(seed,hasher(x));
            return seed;
        }
    };
    template<>
    struct hash<ln::arcID_type>
    {
        // hash function needed to store arcID_type into std::unordered_map
        
        std::size_t operator()(const ln::arcID_type& id)const
        {
            return ln::hash_combine(std::hash<std::uint64_t>()(id.short_channel_id),
                                std::hash<std::uint8_t>()(id.part));
        }
    };
}

namespace ln
{
    // default graph type used for Pickhardtpayments 
    typedef digraph<nodeID_type,arcID_type> graph_type;
    
    // default node and arc index type
    typedef graph_type::node_pos_t node_index_type;
    typedef graph_type::arc_pos_t  arc_index_type;
    
    // default data type to represent costs, capacities and flows
    typedef std::int64_t costflow_value_type;
    
    class base_MCFNetwork
    {
        protected:
        

        
        graph_type G;
        
        std::vector<costflow_value_type> default_cost;
        
        std::vector<costflow_value_type> my_capacity,
                                         my_cost,
                                         my_residual_capacity;
        
        std::vector<costflow_value_type> my_balance;
        
        std::vector<arc_index_type> my_Sarcs,
                                    my_Tarcs;
        
        const node_index_type Source, Sink;
        
        
        void check_data_vectors_size()
        {
            my_capacity.resize(G.max_num_arcs());
            my_cost.resize(G.max_num_arcs());
            default_cost.resize(G.max_num_arcs());
            my_residual_capacity.resize(G.max_num_arcs());
            
            my_balance.resize(G.max_num_nodes());
            
            my_Sarcs.resize(G.max_num_nodes());
            my_Tarcs.resize(G.max_num_nodes());
        }
        
        auto cost_at(arc_index_type arc)const
        {
            return my_cost.at(arc);
        }
        auto flow_at(arc_index_type arc)const
        {
            auto dual = G.arc_dual(arc);
            return my_residual_capacity.at(dual);
        }
        auto capacity_at(arc_index_type arc)const
        {
            return my_capacity.at(arc);
        }
        auto residual_capacity_at(arc_index_type arc)const
        {
            return my_residual_capacity.at(arc);
        }
   
        bool source_is_balanced()const
        {
            bool result=true;
            for(auto arc : G.out_arcs(Source))
                result = result && residual_capacity_at(arc)==0;
            return result;
        }
        bool sink_is_balanced()const
        {
            bool result=true;
            for(auto arc : G.out_arcs(Sink))
                result = result && residual_capacity_at(G.arc_dual(arc))==0;
            return result;
        }
        bool all_nodes_are_balanced()const
        {
            bool result=true;
            for(auto node : G.nodes())
            if(node!=Source && node!=Sink)
            {
                costflow_value_type balance=0;
                for(auto arc: G.out_arcs(node))
                {
                    costflow_value_type r = G.is_dual_arc(arc) ?     
                        residual_capacity_at(arc) : 0-residual_capacity_at(G.arc_dual(arc));
                    balance += r;
                }
                
                result = result && balance==0;
            }
            return result;
        }
        bool node_exists(const nodeID_type id)const
        {
            return G.get_node(id) != G.NONE;
        }
        bool arc_exists(const arcID_type id)const
        {
            return G.get_arc(id) != G.NONE;
        }
        node_index_type check_node(const nodeID_type id)
        {
            if(node_exists(id))
                return G.get_node(id);
            
            return AddNode(id);
        }
        
        void Clear()
        {
            std::fill(my_balance.begin(),my_balance.end(),0);
            std::fill(my_residual_capacity.begin(),my_residual_capacity.end(),0);
            for(auto arc : G.arcs())
                if(!G.is_dual_arc(arc))
                    my_residual_capacity.at(arc) = capacity_at(arc);
        }
   
        public:
        
        /* Initialize the network */
        
        base_MCFNetwork():
            Source(G.new_node()),
            Sink(G.new_node())
        {
        }
        
        node_index_type AddNode(const nodeID_type id)
        {   
            if(node_exists(id))
                return G.get_node(id);
            
            const auto index = G.add_node(id);
            const auto [s_arc,s_dual] = G.add_arc(Source,index); 
            const auto [t_arc,t_dual] = G.add_arc(index,Sink);
            
            check_data_vectors_size();
            
            my_cost.at(s_arc) = my_cost.at(s_dual) = my_cost.at(t_arc) = my_cost.at(t_dual) = 0;
            default_cost.at(s_arc) = default_cost.at(s_dual) = default_cost.at(t_arc) = default_cost.at(t_dual) = 0;
            my_capacity.at(s_dual) = my_capacity.at(t_dual)=0;
            my_capacity.at(s_arc) = my_capacity.at(t_arc)=0;
            
            my_Sarcs.at(index) = s_arc;
            my_Tarcs.at(index) = t_arc;
            return index;
        }
        
        void SetNodeSupply(const nodeID_type id, const costflow_value_type supply)
        {
            const auto index = check_node(id);
            
            const auto s_arc =  my_Sarcs.at(index);
            const auto t_arc =  my_Tarcs.at(index);
            const auto s_dual = G.arc_dual(s_arc);
            const auto t_dual = G.arc_dual(t_arc);
            
            my_cost.at(s_arc) = my_cost.at(t_arc) = my_cost.at(s_dual) = my_cost.at(t_dual) = 0;
            default_cost.at(s_arc) = default_cost.at(t_arc) = default_cost.at(s_dual) = default_cost.at(t_dual) = 0;
            my_capacity.at(s_arc) = my_capacity.at(t_arc) = my_capacity.at(s_dual) = my_capacity.at(t_dual) = 0;
            
            if(supply>0)
            {
                my_capacity.at(s_arc) = supply;
            }
            if(supply<0)
            {
                my_capacity.at(t_arc) = -supply;
            }
        }
       
        std::uint64_t AddArc(const nodeID_type tail, 
                    const nodeID_type head, 
                    const arcID_type id,
                    const costflow_value_type capacity,
                    const costflow_value_type cost)
        {
            if(arc_exists(id))
            {
                throw mcf_arc_duplicate();
            }
            check_node(tail);
            check_node(head);
            
            const auto [arc,dual] = G.add_arc(tail, head, id);
            
            check_data_vectors_size();
            
            my_capacity.at(arc) = capacity;
            my_capacity.at(dual) = 0;
            
            my_cost.at(arc)=cost;
            my_cost.at(dual)=-cost;
            
            default_cost.at(arc)=cost;
            default_cost.at(dual)=-cost;
            
            return arc;
        }
       
        /* Compute */
        virtual void Solve() = 0;
        
        /* Extract information */
       
        bool FlowIsFeasible()const
        {
            return source_is_balanced() && sink_is_balanced() && all_nodes_are_balanced();
        }
        
        costflow_value_type Flow(const arcID_type a)const 
        {
            const auto arc = G.get_arc(a);
            if(arc==G.NONE)
                throw mcf_arc_missing();
            const auto dual = G.arc_dual(arc);
            return my_residual_capacity.at(dual);
        }
        
        int CountFlow()const
        {
            int count=0;
            for(const auto arc : G.arcs())
            {
                const auto dual = G.arc_dual(arc);
                const auto [a,b] = G.arc_ends(arc);
                if(a==Source || a==Sink || b==Source || b==Sink)
                    continue;
                
                const costflow_value_type flow = my_residual_capacity.at(dual);
                if(!G.is_dual_arc(arc) && flow>0)
                {
                    ++count;
                }
            }
            return count;
        }
        
        void GetFlow(std::uint64_t * const index_list, costflow_value_type * const flow_list)const
        {
            int count=0;
            for(const auto arc : G.arcs())
            {
                const auto dual = G.arc_dual(arc);
                const auto [a,b] = G.arc_ends(arc);
                if(a==Source || a==Sink || b==Source || b==Sink)
                    continue;
                    
                const costflow_value_type flow = my_residual_capacity.at(dual);
                if(!G.is_dual_arc(arc) && flow>0)
                {
                    index_list[count]=arc;
                    flow_list[count]=flow;
                    ++count;
                }
            }
        }
        
        /* Update the network */
        
        void UpdateArc(
                    const arcID_type id,
                    const costflow_value_type capacity,
                    const costflow_value_type cost)
        {
            if(!arc_exists(id))
            {
                throw mcf_arc_missing();
            }
            const auto arc = G.get_arc(id);
            const auto dual = G.arc_dual(arc);
            
            my_capacity.at(arc) = capacity;
            my_capacity.at(dual) = 0;
            
            my_cost.at(arc)=cost;
            my_cost.at(dual)=-cost;
        }
        
        void RemoveArc(const arcID_type id)
        {
            if(!arc_exists(id))
            {
                throw mcf_arc_missing();
            }
            G.remove_arc(id);
        }
        
        void forget()
        {
            std::copy(default_cost.begin(),default_cost.end(),my_cost.begin());
            
            for(const auto index : G.nodes())
            if(index!=Source && index!=Sink)
            {
                const auto s_arc =  my_Sarcs.at(index);
                const auto t_arc =  my_Tarcs.at(index);
                const auto s_dual = G.arc_dual(s_arc);
                const auto t_dual = G.arc_dual(t_arc);
            
                my_capacity.at(s_arc) = my_capacity.at(t_arc) = my_capacity.at(s_dual) = my_capacity.at(t_dual) = 0;
            }
        }
        
        virtual ~base_MCFNetwork(){}
    };
    
    class MCFNetwork : public base_MCFNetwork
    {
        using base_MCFNetwork::Clear;
        using base_MCFNetwork::FlowIsFeasible;
        
        // template specialization for Edmond-Karp's algorithm (successive shortest path)
        using MCFAugmentingPaths =
            ln::mincostflow_EdmondsKarp< costflow_value_type, 
                                     ln::shortestPath_FIFO<costflow_value_type> > ;
        
        // template specialization for Goldberg-Tarjan's algorithm (cost-scaling)
        // this template takes a maxflow function as argument, but whatever maxflow we use 
        // it will barely affect the runtime.
        using MCFCostScaling =
            ln::mincostflow_costScaling< 
                ln::maxflow_augmenting_path< costflow_value_type,
                                             ln::pathSearch_labeling> >;
        
        public:
        
        virtual void Solve() override
        // The default solver uses MCFAugmentingPaths
        {
            Solve_by_AugmentingPaths();
        }

        void Solve_by_AugmentingPaths()
        {
            Clear();
            MCFAugmentingPaths f;
            f.solve(G,Source,Sink,my_cost,my_residual_capacity);
            if(!FlowIsFeasible())
                throw mcf_unfeasible();
        }
        void Solve_by_CostScaling()
        // This is a test solver that uses a static graph but with very efficient adjacency access,
        // I put it here because it is faster than MCFCostScaling.
        // Otherwise this fuction could have been written simpy as
        // {
        //     Clear();
        //     MCFCostScaling f;
        //     f.solve(G,Source,Sink,my_cost,my_residual_capacity);
        //     if(!FlowIsFeasible())
        //         throw mcf_unfeasible();
        // }
        {   
            Clear();
            
            const int N = G.num_nodes();
            // const int M = G.num_arcs()/2;
            
            std::vector<costflow_value_type> capacity;
            std::vector<costflow_value_type> cost;
            
            std::vector<int> nodeindex_to_i(G.max_num_nodes());
            // std::vector<node_index_type> i_to_nodeindex(N);
            
            std::vector< std::pair<int,int> > edges;
           
            unsigned int ind=0;
            for(const auto node : G.nodes())
            {
                nodeindex_to_i.at(node)=ind;
                // i_to_nodeindex[i]=node;
                ++ind;
            }
            
            // std::vector<int> arcindex_to_i(G.max_num_arcs());
            std::vector<arc_index_type> i_to_arcindex;
            
            // ind=0;
            for(const auto arc : G.arcs())
            if(!G.is_dual_arc(arc))
            {
                // arcindex_to_i[arc]=i;
                i_to_arcindex.push_back(arc);
                
                cost.push_back ( my_cost.at(arc) );
                capacity.push_back( my_capacity.at(arc) );
                
                const auto [a,b] = G.arc_ends(arc);
                edges.push_back({nodeindex_to_i.at(a),nodeindex_to_i.at(b)});
                
                // ++ind;
            }
            
            simple_mcf<costflow_value_type> f(N,edges,capacity,cost);
            
            auto flow = f.solve(nodeindex_to_i.at(Source),nodeindex_to_i.at(Sink));
            
            assert(edges.size()==flow.size());
            
            for(unsigned int i=0;i<flow.size();++i)
            {
                const auto arc = i_to_arcindex.at(i);
                const auto dual = G.arc_dual(arc);
                my_residual_capacity.at( dual ) = flow.at(i);
                my_residual_capacity.at( arc ) = my_capacity.at(arc) - flow.at(i);
            }
            if(!FlowIsFeasible())
                throw mcf_unfeasible(); 
        }
    };
}

