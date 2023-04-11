
## Fundamental Types

```
typedef int64_t costflow_value_type;
```
The value of costs and capacities are represented as 64 bit integers.

## Node ID

```
struct nodeID_type
{
    uint8_t k[33];
};
```

The node id here is represented as a 33 byte array representing its public key.

## Arc ID

```
struct arcID_type
{
    uint64_t short_channel_id;
    uint8_t part; // 1 bit for the direction 7 bits for the linearization part
};
```

This struct encodes the arc ids with two variables:

```
short_channel_id (8 bytes) = block_height (3 bytes) || transaction_index (3 bytes) || output_index (2 bytes)
```

```
part (1 byte) = piece_id (7 bits) || direction (1 bit) 
```

The `piece_id` is an ordinal number to identify each of the multiple pieces in which a channel is
split in the piecewise linear approximation of the cost function.
The `direction` bit purpose is to distinguish the channel direction.
For a channel that connects nodes `A` and `B`, we have direction equals 0 if `A < B`.
and direction equals 1 if `A>B`.

## MCFNetwork

The struct `MCFNetwork` encondes the state of the graph internally and access the methods for
solving the Minimum Cost Flow problem.

```
MCFNetwork* MCFNetwork_new();
```

Creates a new `MCFNetwork` object.

```
int MCFNetwork_SetNodeSupply(
    MCFNetwork* g,
    const nodeID_type id, 
    const costflow_value_type supply);
```

Sets a node supply/demand of the node `id` with the value `supply`.
Notice that a necessary condition for a feasible flow
is that the sum of all node balances equals to zero.
A new node will be added if its id is new to the graph.
By default node supply is set to 0.

```
int MCFNetwork_AddArc(
    MCFNetwork* g,
    const nodeID_type tail, 
    const nodeID_type head, 
    const arcID_type id,
    const costflow_value_type capacity,
    const costflow_value_type cost,
    uint64_t * p_arc_index);
```

Adds a new arc from nodes tail to head.
The nodes are created if they don't exist.
This function would not allow the creation of existing arcs.

```
int MCFNetwork_Solve(MCFNetwork* g);
```

Solve the Min. Cost Flow problem with a function call. The solution is kept internally.
It will throw an exception if no feasible flow is possible for the given network.


```
int MCFNetwork_Flow(
    MCFNetwork const * g,
    const arcID_type a,
    costflow_value_type* p_value);
```

Returns the value of the flow given to the arc defined by the pair `(arc_id,arc_part)`.
One would call this function after `Solve()`.
It will throw an exception if the arc does not exist.


```
int MCFNetwork_CountFlow(
    MCFNetwork const* g,
    uint64_t* p_n);
```

Returns the number of arcs with non-zero flow.

```
int MCFNetwork_GetFlow(
    MCFNetwork const * g,
    uint64_t * const index_list, 
    costflow_value_type * const flow_list);
```

Returns a list of the arcs that have non-zero flow
and their respective flow values.
This could be faster than calling `Flow(...)` for every
single arc.
Arcs are identified by their internal indexes.

```
int MCFNetwork_UpdateArc(
    MCFNetwork* g,
    const arcID_type id,
    const costflow_value_type capacity,
    const costflow_value_type cost);
```

Updates the parameters of the arc identified by the pair `(arc_id,arc_part)`.
It will throw an exception if the arc does not exist.


```
int MCFNetwork_RemoveArc(
    MCFNetwork*g, 
    const arcID_type id);
```

Removes the arc identified by the pair `(arc_id,arc_part)`.
It will throw an exception if the arc does not exist.


```
int MCFNetwork_forget(MCFNetwork *g);
```

It sets the cost of every arc to their initial value, ie. the value set by `AddArc` 
and it sets the supply/demand of every node to zero.

```    
void MCFNetwork_free(MCFNetwork* g);
```

It frees the allocated memory for the `MCFNetwork` object.

