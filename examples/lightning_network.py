import MinCostFlow

RENE = "03efccf2c383d7bf340da9a3f02e2c23104a0e4fe8ac1a880c8e2dc92fbdacd9df"
C_OTTO = "027ce055380348d7812d2ae7745701c9f93e70c1adeb2657f053f91df4f2843c71"
n1 = "03434a39cd9a537c852fc8fb72454086d726f9111e9f730cef4985c39c11fae944"
n2 =  "027ccec61f4bf1fafb5156931da6527dc104ec3613dd4f4050161d89dd76ab494c"
n3 = "026726a4b043d413b45b334876d17b8a98848129604429ec65532ba286a42efeac"

nodes = [ RENE, C_OTTO, n1, n2 ,n3]

arcs = [ {"id": "1x0x0", "part": 0, "capacity": 2, "cost": 1, "source": RENE, "destination": n1}, 
         {"id": "2x0x0", "part": 0, "capacity": 2, "cost": 3, "source": RENE, "destination": n2},
         {"id": "3x0x0", "part": 0, "capacity": 2, "cost": 5, "source": n1,   "destination": n3},
         {"id": "4x0x0", "part": 0, "capacity": 2, "cost": 1, "source": n2,   "destination": n3},
         {"id": "5x0x0", "part": 0, "capacity": 2, "cost": 0, "source": n3,   "destination": C_OTTO},
         ] 

G = MinCostFlow.MCFNetwork()

for n in nodes:
    G.SetNodeSupply(n,0)
G.SetNodeSupply(C_OTTO,-2)
G.SetNodeSupply(RENE,2)

for arc in arcs:
    G.AddArc(arc['source'],arc['destination'],arc['id'],arc['part'],arc['capacity'],arc['cost'])

G.Solve()


for arc in arcs:
    f= G.Flow(arc['id'],arc['part'])
    print(arc['id'],arc['part'],f)
