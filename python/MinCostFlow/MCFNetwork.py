from .MCF_common import *

## Initialize the Network 

libmcf.MCFNetwork_new.restype = ctypes.c_void_p
# libmcf.MCFNetwork_new.argtypes = 

libmcf.MCFNetwork_SetNodeSupply.restype = ctypes.c_int
libmcf.MCFNetwork_SetNodeSupply.argtypes = [ctypes.c_void_p,
                                            nodeID_type,
                                            ctypes.c_longlong]

libmcf.MCFNetwork_AddArc.restype = ctypes.c_int
libmcf.MCFNetwork_AddArc.argtypes = [ctypes.c_void_p,
                                     nodeID_type,
                                     nodeID_type,
                                     arcID_type,
                                     ctypes.c_longlong, 
                                     ctypes.c_longlong,
                                     ctypes.POINTER(ctypes.c_longlong)]

## Compute

libmcf.MCFNetwork_Solve.restype = ctypes.c_int
libmcf.MCFNetwork_Solve.argtypes =[ctypes.c_void_p]

# use only for testing
libmcf.MCFNetwork_Solve_by_AugmentingPaths.restype = ctypes.c_int
libmcf.MCFNetwork_Solve_by_AugmentingPaths.argtypes =[ctypes.c_void_p]

# use only for testing
libmcf.MCFNetwork_Solve_by_CostScaling.restype = ctypes.c_int
libmcf.MCFNetwork_Solve_by_CostScaling.argtypes =[ctypes.c_void_p]

## Extract information

libmcf.MCFNetwork_Flow.restype = ctypes.c_int
libmcf.MCFNetwork_Flow.argtypes = [ ctypes.c_void_p,
                                    arcID_type,
                                    ctypes.POINTER(ctypes.c_longlong)]

libmcf.MCFNetwork_CountFlow.restype = ctypes.c_int
libmcf.MCFNetwork_CountFlow.argtypes = [ ctypes.c_void_p,
                                         ctypes.POINTER(ctypes.c_longlong)]

libmcf.MCFNetwork_GetFlow.restype= ctypes.c_int
libmcf.MCFNetwork_GetFlow.argtypes = [ctypes.c_void_p,
                                      ctypes.POINTER(ctypes.c_longlong),
                                      ctypes.POINTER(ctypes.c_longlong)]

## Update the network

libmcf.MCFNetwork_UpdateArc.restype= ctypes.c_int
libmcf.MCFNetwork_UpdateArc.argtypes = [ctypes.c_void_p,
                                        arcID_type,
                                        ctypes.c_longlong, 
                                        ctypes.c_longlong]

libmcf.MCFNetwork_RemoveArc.restypes= ctypes.c_int
libmcf.MCFNetwork_RemoveArc.argtypes = [ctypes.c_void_p,
                                        arcID_type]

libmcf.MCFNetwork_forget.restype= ctypes.c_int
libmcf.MCFNetwork_forget.argtypes = [ctypes.c_void_p ]


## Cleanup

# libmcf.MCFNetwork_free.restype
libmcf.MCFNetwork_free.argtypes=[ctypes.c_void_p]

## TODO: proper python error handling

class MCFNetwork(object):   
    '''
    This class is used to represent the Network, which is a directed graph with capacities and costs
    associated to arcs.
    Nodes and arcs are referenced by their respective ids.
    '''
    ## Initialize the network
    def __init__(self):
        self.obj = ctypes.c_void_p(libmcf.MCFNetwork_new())
        self.error_code = 0
    
    def SetNodeSupply(self, node: str, balance: int):
        '''
        Sets a node supply/demand aka balance.
        The convention is that it is a supply if balance>0
        and a demand if balance<0. 
        
        Notice that a necessary condition for a feasible flow
        is that the sum of all node balances = 0.
        
        A new node will be added if its id is new to the graph.
        
        By default node supply is set to 0.
        '''
        self.error_code = libmcf.MCFNetwork_SetNodeSupply(self.obj,
            nodeID_type.from_str(node),
            ctypes.c_longlong(balance))
        
        if self.error_code!=0:
            raise Exception()
    
    def AddArc(self,tail: str, head: str, arc_id: str, arc_part: int, capacity: int, cost: int):
        '''
        Add a new arc from nodes tail to head.
        The nodes are created if they don't exist.
        
        This function would not allow the creation of existing arcs.
        '''
        aid = arcID_type.from_args(arc_id,arc_part)
        
        arc_index = ctypes.c_longlong(0)
        ptr_arc_index = ctypes.pointer(arc_index)
        
        self.error_code = libmcf.MCFNetwork_AddArc(self.obj,
            nodeID_type.from_str(tail),
            nodeID_type.from_str(head),
            aid,
            ctypes.c_longlong(capacity),
            ctypes.c_longlong(cost),
            ptr_arc_index)
        
        if self.error_code!=0:
            raise Exception()
        
        return arc_index.value
   
    ## Compute
    
    def Solve(self):
        '''
        Solve the Min. Cost Flow problem.
        It will throw an exception if no feasible flow is possible for the given network.
        '''
        self.error_code = libmcf.MCFNetwork_Solve(self.obj)
        if self.error_code!=0:
            raise Exception()
    
    def _Solve_by_AugmentingPaths(self):
        '''
        Solve the Min. Cost Flow problem using the Augmenting-paths algorithm.
        It will throw an exception if no feasible flow is possible for the given network.
        Use only for benchmarks.
        '''
        self.error_code = libmcf.MCFNetwork_Solve_by_AugmentingPaths(self.obj)
        if self.error_code!=0:
            raise Exception()
    
    def _Solve_by_CostScaling(self):
        '''
        Solve the Min. Cost Flow problem using the Cost-scaling algorithm.
        It will throw an exception if no feasible flow is possible for the given network.
        Use only for benchmarks.
        '''
        self.error_code = libmcf.MCFNetwork_Solve_by_CostScaling(self.obj)
        if self.error_code!=0:
            raise Exception()
  
    ## Extract information
    
    def Flow(self,arc_id: str, arc_part: int):
        '''
        Returns the value of the flow given to the arc defined by the pair '(arc_id,arc_part)'.
        One would call this function after 'Solve()'
        It will throw an exception if the arc does not exist.
        '''
        aid = arcID_type.from_args(arc_id,arc_part)
        value = ctypes.c_longlong(0)
        ptr_value = ctypes.pointer(value)
        self.error_code = libmcf.MCFNetwork_Flow(self.obj,aid,ptr_value)
        if self.error_code!=0:
            raise Exception()
        return value.value
        
    
    def FlowArray(self):
        '''
        Returns a list of the arcs that have non-zero flow
        and their respective flow values.
        This could be faster than calling `Flow(...)` for every
        single arc.
        Arcs are identified by their internal indexes.
        '''
        N = ctypes.c_longlong(0)
        ptr_N = ctypes.pointer(N)
        
        self.error_code = libmcf.MCFNetwork_CountFlow(self.obj,ptr_N)
        
        if self.error_code!=0:
            raise Exception()
        
        index_list = numpy.arange(N.value, dtype=ctypes.c_longlong)
        flow_list  = numpy.arange(N.value, dtype=ctypes.c_longlong)
        
        self.error_code = libmcf.MCFNetwork_GetFlow(self.obj,
                                  index_list.ctypes.data_as(ctypes.POINTER(ctypes.c_longlong)),
                                  flow_list.ctypes.data_as(ctypes.POINTER(ctypes.c_longlong)),)
        
        if self.error_code!=0:
            raise Exception()
        
        return index_list,flow_list
    
    ## Update the network
    
    def RemoveArc(self,arc_id: str, arc_part: int):
        '''
        Removes the arc identified by the pair '(arc_id,arc_part)'.
        It will throw an exception if the arc does not exist.
        '''
        aid = arcID_type.from_args(arc_id,arc_part)
        self.error_code = libmcf.MCFNetwork_RemoveArc(self.obj,aid)
        
        if self.error_code!=0:
            raise Exception()
        
    def UpdateArc(self,arc_id: str, arc_part: int, capacity: int, cost: int):
        '''
        Removes the parameters of the arc identified by the pair '(arc_id,arc_part)'.
        It will throw an exception if the arc does not exist.
        '''
        aid = arcID_type.from_args(arc_id,arc_part)
        self.error_code = libmcf.MCFNetwork_UpdateArc(self.obj,
            aid,
            ctypes.c_longlong(capacity),
            ctypes.c_longlong(cost))
        
        if self.error_code!=0:
            raise Exception()
    
    def Forget(self):
        '''
        It sets the cost of every arc to their initial value, ie. the value set by `AddArc` 
        and it sets the supply/demand of every node to zero.
        '''
        self.error_code = libmcf.MCFNetwork_forget(self.obj)
        
        if self.error_code!=0:
            raise Exception()
    
    ## Cleanup
   
    def __del__(self):
        libmcf.MCFNetwork_free(self.obj)
    
