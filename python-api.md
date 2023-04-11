# Python API

The library provides a python class `MCFNetwork` which represents the cost network and it contains
methods to compute an optimal solution.

### Arc ID

```
class arcID_type(ctypes.Structure):
```

This class encodes the arc ids with two variables:

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


```
def from_args(arc_id: str, arc_part: int):
```

A function to construct an `arcID_type` from a
short channel id and `arc_part` index.

```
def serialize(self):
```

Returns the short channel id (as a string) and the index of the channel part of the
piecewise linear-cost channels decomposition

### Node ID

```
class nodeID_type(ctypes.Structure):
```

The node id here is represented as a 33 byte array 
representing its public key.

```
def from_str(nodeid : str):
```

It returns a `nodeID_type` from the string that represents the node public key.

```
def serialize(self):
```

It returns the string that encodes the node public key.


### MCFNetwork

```
class MCFNetwork(object):   
```

This class is used to represent the Network, which is a directed graph with capacities and costs
associated to arcs.
Nodes and arcs are referenced by their respective ids.


```
def SetNodeSupply(self, node: str, balance: int):
```

Sets a node supply/demand of the node `node` with the value `balance`.
The convention is that it is a supply if `balance`>0 and a demand if `balance`<0. 
Notice that a necessary condition for a feasible flow
is that the sum of all node balances equals to zero.
A new node will be added if its id is new to the graph.
By default node supply is set to 0.

```
def AddArc(self,tail: str, head: str, arc_id: str, arc_part: int, capacity: int, cost: int):
```
        
Adds a new arc from nodes tail to head.
The nodes are created if they don't exist.
This function would not allow the creation of existing arcs.
    
```    
def Solve(self):
```

Solve the Min. Cost Flow problem with a function call. The solution is kept internally.
It will throw an exception if no feasible flow is possible for the given network.

```     
def Flow(self,arc_id: str, arc_part: int):
```

Returns the value of the flow given to the arc defined by the pair `(arc_id,arc_part)`.
One would call this function after `Solve()`.
It will throw an exception if the arc does not exist.
        
```    
def FlowArray(self):
```

Returns a list of the arcs that have non-zero flow
and their respective flow values.
This could be faster than calling `Flow(...)` for every
single arc.
Arcs are identified by their internal indexes.
    
```
def RemoveArc(self,arc_id: str, arc_part: int):
```

Removes the arc identified by the pair `(arc_id,arc_part)`.
It will throw an exception if the arc does not exist.
        
```
def UpdateArc(self,arc_id: str, arc_part: int, capacity: int, cost: int):
```
        
Updates the parameters of the arc identified by the pair `(arc_id,arc_part)`.
It will throw an exception if the arc does not exist.
        
```    
def Forget(self):
```
It sets the cost of every arc to their initial value, ie. the value set by `AddArc` 
and it sets the supply/demand of every node to zero.
