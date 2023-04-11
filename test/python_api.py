from MinCostFlow import MCFNetwork
import numpy

my_mcf = MCFNetwork()

KEY_BYTES = 33
node=[ '00' * KEY_BYTES ,
       '11' * KEY_BYTES ,
       '22' * KEY_BYTES ,
       '33' * KEY_BYTES ]

result_flow = numpy.array([1,1,0,0,1,2])
edges = [ (node[0],node[3]), 
          (node[0],node[2]), 
          (node[1],node[2]), 
          (node[1],node[0]),
          (node[2],node[3]),
          (node[3],node[1])]
capacity = [2,1,1,1,4,2]
unit_cost = [4,1,0,1,2,0]
short_channel_id = ['0x0x0','1x1x1','2x2x2','3x3x3','4x4x4','5x5x5']

def test_case(solver,arc_list,short_channel_id_list,cap_list,cost_list,fiducial_result):
    solver.SetNodeSupply(node[0],2)    
    solver.SetNodeSupply(node[1],-2)    
    map_arc={}
    for i in range(len(arc_list)):
        (a,b) = arc_list[i]
        j = solver.AddArc(a,b,short_channel_id_list[i],0,cap_list[i],cost_list[i])
        map_arc[j]=i
    solver.Solve()
    index,sol = solver.FlowArray() 
    flow_sol = numpy.zeros(len(arc_list),dtype=int)
    for j,f in zip(index,sol):
        i = map_arc[j]
        flow_sol[i]=f
    diff = numpy.abs(flow_sol - fiducial_result).max()
    if diff!=0:
        return False
    return True
    
if not test_case(my_mcf,edges,short_channel_id,capacity,unit_cost,result_flow):
    raise Exception('test case failed')
