"""
CVRP using MiniZinc
"""
import sys
from minizinc import Instance, Model, Solver
import numpy as np
from scipy.spatial.distance import cityblock
import datetime
import networkx as nx
import matplotlib.pyplot as plt
import random

# Here I plot the solution if any
def result_printer(result,num_courier,num_item,xc,yc):
    if result.status.has_solution():
        out = open('../out/out'+sys.argv[1][1:][3:]+'.txt', 'w')
        print("Solution found")
        out.write("Solution found")
        out.write("\n")
        print()
        lines = str(result.solution).split("\n")
        for i in range(len(lines)):
            lines[i] = lines[i].strip("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ[]= ;\n").replace(",","")
        objective = int(lines[0])
        print("Total length of the routes: "+ str(objective))
        out.write("Total length of the routes: "+ str(objective))
        out.write("\n")
        vehicle = list(map(int, lines[1].split()))
        successor = list(map(int, lines[2].split()))
        predecessor = list(map(int, lines[3].split()))
        load = list(map(int, lines[4].split()))
        tours = [[] for x in range(1, num_courier+1)]
        loads = [[] for x in range(1, num_courier+1)]
        active_arcs = [[0] for x in range(1, num_courier+1)]
        for c in range(num_courier):
            for i in range(len(vehicle)):
                if vehicle[i] == c+1 and load[i]!=0 and predecessor[i] < num_item+1 :
                    tours[c].append(predecessor[i]-1)
                    loads[c].append(load[i])
        for c in range(num_courier):
            zipped_pairs = zip(tours[c],loads[c])
            zipped_pairs = sorted(zipped_pairs, key = lambda x: x[1])
            active_arcs[c] = [item[0] for item in zipped_pairs]
            active_arcs[c].append(num_item)
            active_arcs[c].insert(0,num_item)
        for c in range(num_courier):
            print("Route courier "+str(c+1)+":")
            out.write("Route courier "+str(c+1)+":")
            out.write("\n")
            tmp="0->"
            for i in range(1,len(active_arcs[c])-1):
                tmp += str(active_arcs[c][i]+1)+'->'
            tmp += str(0)
            print(tmp)
            out.write(tmp)
            out.write("\n")
        out.close()
        
        # Plot the tour of each vehicle
        fig, ax = plt.subplots()
        G=nx.Graph(name="route_mnz")
        edges = []
        for r in active_arcs:
            route_edges = [(r[n],r[n+1]) for n in range(len(r)-1)]
            G.add_nodes_from(r)
            G.add_edges_from(route_edges)
            edges.append(route_edges)
        coordinates = np.array([ [x,y] for x,y in zip(xc, yc)])
        nodes = np.arange(0,num_item+1)
        npos = dict(zip(nodes, coordinates))
        nx.draw_networkx_nodes(G,pos=npos, node_size=9)
        colors = [ "#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)]) for i in range(len(active_arcs))]
        linewidth = 2
        plt.plot(xc[-1], yc[-1], c='r', marker='s')
        for ctr, edgelist in enumerate(edges):
            nx.draw_networkx_edges(G,pos=npos,edgelist=edgelist,edge_color = colors[ctr], width=linewidth)
        plt.axis('on')
        ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=True)
        plt.savefig('../plots/route_mnz_'+sys.argv[1]+'.png')
        return
        
    else:
        print("No solution found before the timeout.")

def main(argv):
    # Create a MiniZinc model
    gecode = Solver.lookup("gecode")

    model=Model()
    model.add_file("HamiltonianCP.mzn")

    # Transform Model into a instance
    inst = Instance(gecode, model)

    # Instantiate variables from file
    data_file= open('./MCP_Instances/'+argv)
    lines=[]
    for line in data_file:
      lines.append(line)
    data_file.close()
    num_courier=int(lines[0].rstrip('\n'))
    num_item=int(lines[1].rstrip('\n'))
    load=list(map(int, lines[2].rstrip('\n').split()))
    weight=list(map(int, lines[3].rstrip('\n').split()))
    x=list(map(int, lines[4].rstrip('\n').split()))
    y=list(map(int, lines[5].rstrip('\n').split()))
    dist=np.zeros((len(x),len(x)))
    for i in range(len(x)):
      for j in range(len(x)):
        dist[i,j] = int(cityblock([x[i],y[i]],[x[j],y[j]]))
    dist=dist.astype(int)

    # Reduction of number couriers
    load = sorted(load, reverse=True)
    for i in range(len(load)):
        if sum(load[0:i])>= sum(weight):
            load1 = load[0:i]
            load=load1
            num_courier = i
            break   
    inst["courier"] = num_courier
    inst["items"] = num_item
    inst["load"] = load
    inst["weight"] = weight
    inst["dist"]= dist

    # Output
    result = inst.solve(timeout=datetime.timedelta(minutes=5))
    result_printer(result, num_courier, num_item, x, y)

if __name__ == "__main__":
    main(sys.argv[1])