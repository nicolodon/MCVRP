"""
CVRP using Z3 and SMT
"""
import sys
import numpy as np
from scipy.spatial.distance import cityblock
import networkx as nx
import matplotlib.pyplot as plt
import random
from z3 import Int,Bool,Solver,Optimize,PbLe,PbEq,Implies,And,Not,Or,If,Sum,sat,set_option

# Save the solution if any   
def on_model(model):
    out = open('../out/out'+sys.argv[1][1:][3:]+'.txt', 'w')
    print("Solution found")
    out.write("Solution found")
    out.write("\n")
    tot = Int('tot')
    print("Total length of the routes: "+ str(model[tot]))
    out.write("Total length of the routes: "+ str(model[tot]))
    out.write("\n")
    track = [[[Bool(f"x_{i}_{j}_{k}") for k in range(num_courier)] for j in range(num_item+1)] for i in range(num_item+1)]
    t = [[[ model.evaluate(track[i][j][k]) for k in range(num_courier)]for j in range(num_item+1)]for i in range(num_item+1)]
    edges=[[] for c in range(num_courier)]
    for c in range(num_courier):
        print("Route courier "+str(c+1)+":")
        out.write("Route courier "+str(c+1)+":")
        out.write("\n")
        active_arcs={}
        nodes=list(range(num_item))
        nodes.insert(0,num_item)
        for i in nodes:
            for j in nodes:
                if t[i][j][c]:
                    active_arcs[i]=j
        finish = False
        i = num_item
        route=str(0)+'->'
        while not finish :
            j=active_arcs.get(i)
            edges[c].append((i,j))
            i=j
            if(j==num_item):
                route+=str(0)
                finish=True
            else:
                route+= str(j+1)+'->'
        print(route)
        out.write(route)
        out.write("\n")
    print()
    out.close()

    # Plot the tour of each vehicle
    fig, ax = plt.subplots()
    G=nx.Graph(name="route_smt")
    for c in range(num_courier):
        G.add_nodes_from(edges[c][:][0])
        G.add_edges_from(edges[c])
    coordinates = np.array([ [x,y] for x,y in zip(xc, yc)])
    nodes = np.arange(0,num_item+1)
    npos = dict(zip(nodes, coordinates))
    nx.draw_networkx_nodes(G,pos=npos, node_size=9)
    colors = [ "#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)]) for i in range(len(edges))]
    linewidth = 2
    plt.plot(xc[-1], yc[-1], c='r', marker='s')
    for ctr, edgelist in enumerate(edges):
        nx.draw_networkx_edges(G,pos=npos,edgelist=edgelist,edge_color = colors[ctr], width=linewidth)
    plt.axis('on')
    ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=True)
    plt.savefig('../plots/route_smt_'+sys.argv[1]+'.png')
  

def courier_scheduling_problem_MTZ(num_items, num_couriers, items_weights, courier_loads, dist):
    track = [[[Bool(f"x_{i}_{j}_{k}") for k in range(num_couriers)] for j in range(num_items+1)] for i in range(num_items+1)]
    u = [Int(f"u_{j}")for j in range(num_items)]
    tot = Int("tot")
    s = Optimize()
    s.set_on_model(on_model)
    s.set("maxsat_engine",'core_maxsat')
    s.set("timeout",300000)
    #s = Solver()

    # Main diagonal equal to 0 
    for i in range(num_items+1):
      #s.add(And([Not(track[i][i][c]) for c in range(num_couriers)]))
      s.add(PbEq([
          (track[i][i][c], 1) for c in range(num_couriers)
        ], 0))
      

    # Each node visited only once
    for j in range(num_items):
        #s.add(And(exactly_one([track[i][j][c] for c in range(num_courier) for i in range(num_items+1)]), exactly_one([track[j][i][c] for c in range(num_courier) for i in range(num_items+1)])))
        s.add(And(PbEq([
                    (track[i][j][c],1) for c in range(num_couriers) for i in range(num_items+1)
                    ],1),
                  PbEq([
                    (track[j][i][c],1) for c in range(num_couriers) for i in range(num_items+1)
                    ],1)))

    
    # Each courier can go back and depart from depot at maximum one time
    for c in range(num_couriers):
        #s.add(And(exactly_one([track[num_items][j][c] for j in range(num_items)]), exactly_one([track[j][num_items][c] for j in range(num_items)])))
        s.add(And(PbEq([(track[num_items][j][c],1) for j in range(num_items)],1),
                  PbEq([(track[j][num_items][c],1) for j in range(num_items)],1)))


    # Constraint on the load of each courier
    for c in range(num_couriers):
      #s.add(courier_loads[c]>= Sum([If(track[i][j][c],items_weights[i],0) for i in range(num_items) for j in range(num_items+1)]))
      s.add(PbLe([
            (track[i][j][c], items_weights[i]) for i in range(num_items) for j in range(num_items+1)
        ], courier_loads[c]))

    # N arcs in = n arcs out
    for c in range(num_couriers):
      for j in range(num_items+1):
        # orig
        s.add(Sum([If(track[i][j][c],1,0) for i in range(num_items+1)])==Sum([If(track[j][i][c],1,0) for i in range(num_items+1)]))
        #s.add(PbEq([(track[i][j][c],1) for i in range(num_items+1)],Sum([(track[j][i][c],1) for i in range(num_items+1)])))


    #  Miller-Tucker-Zemlin formulation
    for c in range(num_courier):
      for i in range(num_items):
        for j in range(num_items):
          s.add(u[i] + If(track[i][j][c],1,0) <= u[j] + num_items*(1-If(track[i][j][c],1,0)))
          s.add(u[i] > 0)
            

    # Constraint to obtain the total length of the routes 
    s.add(tot==Sum([If(track[i][j][c],dist[i][j],0) for i in range(num_items+1) for j in range(num_items+1) for c in range(num_couriers)]))
    s.minimize(tot)

    print("Model loaded, starting to solve the "+sys.argv[1]+".")

    if s.check() == sat:
      m = s.model()
      t = m.evaluate(tot)
      r = [[[m.evaluate(track[i][j][k]) for k in range(num_couriers)]for j in range(num_items+1)]for i in range(num_items+1)]
      for c in range(num_couriers):
         print('Courier', c)
         nodes=list(range(num_items))
         nodes.insert(0,num_items)
         for i in nodes:
           for j in nodes:
              if r[i][j][c] == True:
                print(i,'->',j)
      print("satttttttttttttttttttttttttt")
      print(r)
      return t
    else:
      print("unsat")



def main(argv):
    # Instantiate variables from file
    file= open('./MCP_Instances/'+argv)
    lines=[]
    for line in file:
      lines.append(line)
    file.close()
    global num_courier
    num_courier=int(lines[0].rstrip('\n'))
    global num_item
    num_item=int(lines[1].rstrip('\n'))
    load=list(map(int, lines[2].rstrip('\n').split()))
    weight=list(map(int, lines[3].rstrip('\n').split()))
    global xc
    xc=list(map(int, lines[4].rstrip('\n').split()))
    global yc
    yc=list(map(int, lines[5].rstrip('\n').split()))
    load = sorted(load, reverse=True)
    for i in range(len(load)):
        if sum(load[0:i])>= sum(weight):
            load1 = load[0:i]
            load=load1
            num_courier = i
            break
    dist=np.zeros((len(xc),len(xc)))
    for i in range(len(xc)):
      for j in range(len(xc)):
        dist[i,j] = int(cityblock([xc[i],yc[i]],[xc[j],yc[j]]))
    dist=dist.astype(int)
    dist=dist.tolist()
    # Call to the Z3 SMT MAXSAT solver
    tot = courier_scheduling_problem_MTZ(num_item, num_courier, weight, load, dist)

    """Alternative model that uses lazy constraints: to use it comment the lines above and uncomment the lines below"""
    # tot,best_track = courier_scheduling_problem_MTZ_lazy(num_item, num_courier, weight, load, dist)
    # out = open('../out/out'+sys.argv[1][1:][3:]+'.txt', 'w')
    # print("Solution found")
    # out.write("Solution found")
    # out.write("\n")
    # print("Total length of the routes: "+ str(tot))
    # out.write("Total length of the routes: "+ str(tot))
    # out.write("\n")
    # edges=[[] for c in range(num_courier)]
    # for c in range(num_courier):
    #     print('Courier', c+1)
    #     active_arcs={}
    #     nodes=list(range(num_item))
    #     nodes.insert(0,num_item)
    #     for i in nodes:
    #         for j in nodes:
    #             if best_track[i][j][c]:
    #                 active_arcs[i]=j
    #     finish = False
    #     i = num_item
    #     route=str(i)+'->'
    #     while not finish :
    #         j=active_arcs.get(i)
    #         edges[c].append((i,j))
    #         i=j
    #         if(j==num_item):
    #             route+=str(num_item)
    #             finish=True
    #         else:
    #             route+= str(j)+'->'
    #     print(route)
    #     out.write(route)
    #     out.write("\n")
    # print()
    # out.close()
    # # plot the tour of each vehicle
    # fig, ax = plt.subplots()
    # G=nx.Graph(name="route_smt")
    # for c in range(num_courier):
    #     G.add_nodes_from(edges[c][:][0])
    #     G.add_edges_from(edges[c])
    # coordinates = np.array([ [x,y] for x,y in zip(xc, yc)])
    # nodes = np.arange(0,num_item+1)
    # npos = dict(zip(nodes, coordinates))
    # nx.draw_networkx_nodes(G,pos=npos, node_size=35)
    # colors = [ "#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)]) for i in range(len(edges))]
    # linewidths = 2
    # for ctr, edgelist in enumerate(edges):
    #     nx.draw_networkx_edges(G,pos=npos,edgelist=edgelist,edge_color = colors[ctr], width=linewidth)
    # plt.axis('on')
    # plt.grid(visible=True)
    # ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=True)
    # plt.savefig('../plots/route_smt_lazy_'+sys.argv[1]+'.png')

if __name__ == "__main__":
    main(sys.argv[1])


def courier_scheduling_problem_MTZ_lazy(num_items, num_couriers, items_weights, courier_loads, dist):
    track = [[[Bool(f"x_{i}_{j}_{k}") for k in range(num_courier)] for j in range(num_items+1)] for i in range(num_items+1)]
    u = [Int(f"u_{j}")for j in range(num_items)]
    tot = Int("tot")
    s = Solver()
    s.set("timeout",300000)

    # Tot constraint
    s.add(tot==Sum([If(track[i][j][c],dist[i][j],0) for i in range(num_items+1) for j in range(num_items+1) for c in range(num_couriers)]))
    set_option("verbose", 2)

    # Main diagonal equal to 0 
    for i in range(num_items+1):
      #s.add(And([Not(track[i][i][c]) for c in range(num_couriers)]))
      s.add(PbEq([
          (track[i][i][c], 1) for c in range(num_couriers)
        ], 0))
      

    # Each node visited only once
    for j in range(num_items):
        # orig
        #s.add(And(exactly_one([track[i][j][c] for c in range(num_courier) for i in range(num_items+1)]), exactly_one([track[j][i][c] for c in range(num_courier) for i in range(num_items+1)])))
        s.add(And(PbEq([
                    (track[i][j][c],1) for c in range(num_courier) for i in range(num_items+1)
                    ],1),
                  PbEq([
                    (track[j][i][c],1) for c in range(num_courier) for i in range(num_items+1)
                    ],1)))

    
    # # Each courier can go back and depart from depot at maximum one time
    for c in range(num_courier):
        s.add(And(PbEq([(track[num_items][j][c],1) for j in range(num_items)],1), PbEq([(track[j][num_items][c],1) for j in range(num_items)],1)))

    # Constraint on the load of each courier
    for c in range(num_courier):
      #s.add(courier_loads[c]>= Sum([If(track[i][j][c],items_weights[i],0) for i in range(num_items) for j in range(num_items+1)]))
      s.add(PbLe([
            (track[i][j][c], items_weights[i]) for i in range(num_items) for j in range(num_items+1)
        ], courier_loads[c]))

    # N arcs in = n arcs out
    for c in range(num_courier):
      for j in range(num_items):
        # orig
        s.add(Sum([If(track[i][j][c],1,0) for i in range(num_items+1)])==Sum([If(track[j][i][c],1,0) for i in range(num_items+1)]))
        #s.add(PbEq([(track[i][j][c],1) for i in range(num_items+1)],1))

    # Constraint to eliminate the i->j->i possibility
    for c in range(num_courier):
      for i in range(num_items):
        for j in range(num_items):
          s.add(Implies(track[i][j][c], Not(track[j][i][c])))

    # Constraint to compute the sum of our tour
    best = 1000000
    best_track=[[[False]]]
    for i in range(1000):
      if s.check() == sat:
        m = s.model()
        t = m.evaluate(tot)
        routes = [[ [ m.evaluate(track[i][j][k]) for k in range(num_couriers)] for j in range(num_items+1)] for i in range(num_items+1)]
        subtour=False
        tours =[[] for i in range(num_couriers) ]
        for c in range(num_courier):
          list1 = list(range(num_items))
          list1.insert(0,num_items)
          for i in list1:
            for j in list1:
              if routes[i][j][c] == True:
                tours[c].append((i,j))
        for c in range(num_couriers):
            if(len(tours[c])>0):
              tmp=[tours[c][0]]
              i=1
              while i < len(tours[c]):
                if(tmp[-1][1]==tours[c][i][0]):
                  if tours[c][i] not in tmp:
                    tmp.append(tours[c][i])
                    i=0
                i+=1
              if(len(tmp)!=len(tours[c])):
                subtour=True
                sub=[i[0] for i in [j for j in tmp]]
                for i in sub:
                  for k in range(num_items):
                    if(i!=num_items):
                      s.add(u[i] + If(track[i][k][c],1,0) <= u[k] + num_items*(1-If(track[i][k][c],1,0)))
        if (not subtour) and (t.as_long() < best) :
            s.add(tot < t.as_long())
            best = t.as_long()
            best_track=routes
      else:
        print("unsat")
    return best,best_track
